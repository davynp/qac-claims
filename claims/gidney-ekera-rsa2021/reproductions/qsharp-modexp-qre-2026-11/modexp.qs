// Windowed modular exponentiation of Gidney and Ekerå (2021), Section 2, as a Q# circuit.
//
// Registers use the coset representation of modular integers (Sec. 2.4) and oblivious carry
// runways (Sec. 2.6). A register is stored as pieces: piece i holds its low bits (c_sep bits, fewer
// for the top piece) followed by c_pad padding qubits. For i < w - 1 the padding is the oblivious carry
// runway; for the top piece it is the coset padding. The register's value is
//     sum_i value(piece i) * 2^(i * c_sep),
// and the integer it represents is that value modulo N.
//
// Each windowed multiplication (Sec. 2.5) is two multiply-adds, y += x*k then x += y*(-k^-1), followed
// by a swap of the two registers. A multiply-add iterates over the factor register's bits in windows
// of c_mul bits. For each window it looks up (Babbush et al. QROM, here the QDK's Select) a classically
// precomputed multiple of the constant, addressed by the exponent window and the factor window, adds
// it into the target with a Cuccaro ripple-carry adder confined to each runway piece, and uncomputes
// the lookup with measurement-based uncomputation. Before iterating, each runway is temporarily
// reduced to one carry qubit (Fig. 3, Sec. 2.7). The exponent is consumed by the semi-classical
// Fourier transform with qubit recycling (Sec. 2.3, Fig. 2).

import Std.Arithmetic.*;
import Std.Arrays.*;
import Std.Convert.*;
import Std.Diagnostics.*;
import Std.Math.*;
import Std.Measurement.*;
import Std.ResourceEstimation.*;
import Std.TableLookup.*;

struct Config {
    // Number of bits n of the modulus.
    NumBits : Int,
    // Runway separation c_sep.
    RunwaySep : Int,
    // Padding and runway length c_pad.
    Padding : Int,
    // Multiplication window c_mul.
    MulWindow : Int,
    // Exponent window sizes, in the order the windows are processed (most significant first).
    ExpWindows : Int[],
    // For each exponent window, the classical base b whose power b^e the window multiplies in.
    Bases : BigInt[],
    // The modulus N.
    Modulus : BigInt,
    // Prepare coset padding and runways in superposition (true), or as zero (exact classical check).
    Superpose : Bool,
    // Use all-zero lookup tables (for counting only; the Toffoli count of a lookup does not depend on
    // its contents).
    Placeholder : Bool,
    // Wrap repeated blocks in the QDK's estimate caching (for counting at full size).
    Cache : Bool,
}

function MakeConfig(
    numBits : Int,
    runwaySep : Int,
    padding : Int,
    mulWindow : Int,
    expWindows : Int[],
    bases : BigInt[],
    modulus : BigInt,
    superpose : Bool,
    placeholder : Bool,
    cache : Bool
) : Config {
    new Config {
        NumBits = numBits,
        RunwaySep = runwaySep,
        Padding = padding,
        MulWindow = mulWindow,
        ExpWindows = expWindows,
        Bases = bases,
        Modulus = modulus,
        Superpose = superpose,
        Placeholder = placeholder,
        Cache = cache,
    }
}

// ---------------------------------------------------------------------------------------------
// Register layout
// ---------------------------------------------------------------------------------------------

/// Number of low bits in each piece.
function PieceLowLens(cfg : Config) : Int[] {
    let w = (cfg.NumBits + cfg.RunwaySep - 1) / cfg.RunwaySep;
    mutable lens = [];
    for i in 0..w - 1 {
        set lens += [MinI(cfg.RunwaySep, cfg.NumBits - i * cfg.RunwaySep)];
    }
    lens
}

/// Number of qubits in each piece (low bits plus padding).
function PieceSizes(cfg : Config) : Int[] {
    Mapped(len -> len + cfg.Padding, PieceLowLens(cfg))
}

function SplitRegister(cfg : Config, qs : Qubit[]) : Qubit[][] {
    let sizes = PieceSizes(cfg);
    Partitioned(sizes, qs)[0..Length(sizes) - 1]
}

/// The low bits of every piece, in order of significance.
function LowBits(reg : Qubit[][], lens : Int[]) : Qubit[] {
    mutable bits = [];
    for i in IndexRange(reg) {
        set bits += reg[i][0..lens[i] - 1];
    }
    bits
}

/// The register's bits after its runways have been reduced to single carry qubits (Fig. 3), with the
/// power of two each bit carries. Runway qubits are excluded: their value has been moved into the
/// next piece.
function FactorBits(cfg : Config, reg : Qubit[][], carries : Qubit[]) : (Qubit[], Int[]) {
    let lens = PieceLowLens(cfg);
    let w = Length(reg);
    mutable bits = [];
    mutable weights = [];
    for i in 0..w - 1 {
        for b in 0..lens[i] - 1 {
            set bits += [reg[i][b]];
            set weights += [i * cfg.RunwaySep + b];
        }
        if i < w - 1 {
            set bits += [carries[i]];
            set weights += [(i + 1) * cfg.RunwaySep + cfg.Padding];
        } else {
            for b in 0..cfg.Padding - 1 {
                set bits += [reg[i][lens[i] + b]];
                set weights += [cfg.NumBits + b];
            }
        }
    }
    (bits, weights)
}

// ---------------------------------------------------------------------------------------------
// Cuccaro, Draper, Kutin, and Moulton ripple-carry adder (the paper's Ref. [18])
// ---------------------------------------------------------------------------------------------

operation Maj(c : Qubit, b : Qubit, a : Qubit) : Unit is Adj {
    CNOT(a, b);
    CNOT(a, c);
    CCNOT(c, b, a);
}

operation Uma(c : Qubit, b : Qubit, a : Qubit) : Unit is Adj {
    CCNOT(c, b, a);
    CNOT(a, c);
    CNOT(c, b);
}

/// b += a (mod 2^n). Uses 2(n - 1) Toffolis.
operation CuccaroAdd(a : Qubit[], b : Qubit[]) : Unit is Adj {
    let n = Length(a);
    Fact(n == Length(b), "CuccaroAdd: registers must have equal length");
    if n == 1 {
        CNOT(a[0], b[0]);
    } else {
        use c = Qubit();
        Maj(c, b[0], a[0]);
        for i in 1..n - 2 {
            Maj(a[i - 1], b[i], a[i]);
        }
        CNOT(a[n - 2], b[n - 1]);
        CNOT(a[n - 1], b[n - 1]);
        for i in n - 2..-1..1 {
            Uma(a[i - 1], b[i], a[i]);
        }
        Uma(c, b[0], a[0]);
    }
}

/// b += a, with the carry out XORed into z. Uses 2n Toffolis.
operation CuccaroAddCarry(a : Qubit[], b : Qubit[], z : Qubit) : Unit is Adj {
    let n = Length(a);
    Fact(n == Length(b), "CuccaroAddCarry: registers must have equal length");
    use c = Qubit();
    Maj(c, b[0], a[0]);
    for i in 1..n - 1 {
        Maj(a[i - 1], b[i], a[i]);
    }
    CNOT(a[n - 1], z);
    for i in n - 1..-1..1 {
        Uma(a[i - 1], b[i], a[i]);
    }
    Uma(c, b[0], a[0]);
}

// ---------------------------------------------------------------------------------------------
// Registers: initialization (Sec. 2.15), runway reduction (Fig. 3)
// ---------------------------------------------------------------------------------------------

/// Prepares a register holding `value` in the coset representation with oblivious carry runways.
/// Coset padding: for each padding bit j, a qubit q in |+> controls adding N * 2^j; since the register
/// is then below N * 2^(j+1), q equals [x >= N * 2^j] and is uncomputed by a comparison
/// (ApplyIfLessOrEqualL(X, c, x, q) flips q when c <= x).
/// Runways: each runway is put in uniform superposition r and r is subtracted from the next piece,
/// leaving the represented value unchanged.
operation InitRegister(cfg : Config, reg : Qubit[][], value : BigInt) : Unit {
    let lens = PieceLowLens(cfg);
    let w = Length(reg);
    let flat = LowBits(reg, lens) + reg[w - 1][lens[w - 1]...];
    ApplyXorInPlaceL(value, flat);
    if cfg.Superpose {
        for j in 0..cfg.Padding - 1 {
            use q = Qubit();
            H(q);
            Controlled IncByL([q], (cfg.Modulus <<< j, flat));
            ApplyIfLessOrEqualL(X, cfg.Modulus <<< j, flat, q);  // fires when N * 2^j <= x
        }
        // Top piece first, so that a borrow from a subtraction lands in a runway that is already in
        // superposition.
        for i in w - 2..-1..0 {
            let runway = reg[i][lens[i]...];
            ApplyToEach(H, runway);
            Adjoint IncByLE(runway, reg[i + 1]);
        }
    }
}

/// Fig. 3: add each runway into the low c_pad bits of the next piece, with the carry into a fresh
/// qubit. Afterwards the runway qubits no longer contribute to the represented value.
operation ReduceRunways(cfg : Config, reg : Qubit[][], carries : Qubit[]) : Unit is Adj {
    let lens = PieceLowLens(cfg);
    for i in 0..Length(reg) - 2 {
        CuccaroAddCarry(reg[i][lens[i]...], reg[i + 1][0..cfg.Padding - 1], carries[i]);
    }
}

// ---------------------------------------------------------------------------------------------
// Lookup addition (Sec. 2.5, 2.8)
// ---------------------------------------------------------------------------------------------

/// Lookup table for one window of factor bits: entry (e, v) = k_e * sum_b v_b 2^(weight_b) mod N,
/// indexed by v + 2^(window) * e, as little-endian bits of length n.
function LookupTable(ks : BigInt[], weights : Int[], modulus : BigInt, numBits : Int) : Bool[][] {
    let c = Length(weights);
    mutable data = [];
    for ke in ks {
        for v in 0..(1 <<< c) - 1 {
            mutable s = 0L;
            for b in 0..c - 1 {
                if ((v >>> b) &&& 1) == 1 {
                    set s += 1L <<< weights[b];
                }
            }
            set data += [BigIntAsBoolArray((ke * (s % modulus)) % modulus, numBits)];
        }
    }
    data
}

/// One lookup addition: look up the table entry addressed by (factor window, exponent window), add it
/// into every piece of the target, and uncompute the lookup (measurement-based, inside Select).
operation LookupAdd(data : Bool[][], address : Qubit[], addend : Qubit[][], lens : Int[], target : Qubit[][]) : Unit {
    let out = LowBits(addend, lens);
    Select(data, address, out);
    for i in IndexRange(target) {
        CuccaroAdd(addend[i], target[i]);
    }
    Adjoint Select(data, address, out);
}

/// target += factor * k_e (mod N, approximately), where e is the value of the exponent window.
operation MultiplyAdd(cfg : Config, expWindow : Qubit[], ks : BigInt[], factor : Qubit[][], target : Qubit[][], addend : Qubit[][]) : Unit {
    let lens = PieceLowLens(cfg);
    use carries = Qubit[Length(factor) - 1];
    within {
        ReduceRunways(cfg, factor, carries);
    } apply {
        let (bits, weights) = FactorBits(cfg, factor, carries);
        for start in 0..cfg.MulWindow..Length(bits) - 1 {
            let stop = MinI(start + cfg.MulWindow, Length(bits)) - 1;
            let address = bits[start..stop] + expWindow;
            let tableLen = 1 <<< Length(address);
            if not cfg.Cache or BeginEstimateCaching($"lookup_add_a{Length(address)}", 0) {
                let data = cfg.Placeholder ? [[false, size = cfg.NumBits], size = tableLen] | LookupTable(ks, weights[start..stop], cfg.Modulus, cfg.NumBits);
                LookupAdd(data, address, addend, lens, target);
                if cfg.Cache {
                    EndEstimateCaching();
                }
            }
        }
    }
}

/// Powers base^e mod N for every value e of a window, and the matching -(base^e)^-1 mod N.
function WindowConstants(base : BigInt, windowSize : Int, modulus : BigInt) : (BigInt[], BigInt[]) {
    mutable ks = [];
    mutable kinv = [];
    mutable p = 1L;
    for e in 0..(1 <<< windowSize) - 1 {
        set ks += [p];
        set kinv += [(modulus - InverseModL(p, modulus)) % modulus];
        set p = (p * base) % modulus;
    }
    (ks, kinv)
}

// ---------------------------------------------------------------------------------------------
// The modular exponentiation
// ---------------------------------------------------------------------------------------------

/// Semi-classical inverse Fourier transform on one window (Griffiths and Niu, Sec. 2.3): qubits are
/// measured from the most significant one down; each is rotated by the phase implied by the bits
/// already measured, then measured in the X basis. When counting, the angle depends on measurement
/// outcomes the counter does not know, so one generic rotation per qubit is applied instead (Sec. 2.3:
/// the phase shifts on a qubit combine into a single shift).
operation MeasureWindowSemiclassical(window : Qubit[], measured : Bool[], counting : Bool) : Bool[] {
    mutable ms = measured;
    for q in Reversed(window) {
        let k = Length(ms);
        mutable angle = 0.0;
        for l in MaxI(0, k - 60)..k - 1 {
            if ms[l] {
                set angle += 2.0 * PI() / IntAsDouble(1 <<< (k - l + 1));
            }
        }
        R1(counting ? 0.1 | -angle, q);
        H(q);
        set ms += [MResetZ(q) == One];
    }
    ms
}

/// Runs the modular exponentiation. With a classical exponent (`exponentBits`, least significant bit
/// first), returns the measured accumulator and work registers. With `quantumExponent`, the exponent
/// starts in uniform superposition and is consumed by the semi-classical Fourier transform; the
/// measured frequency bits are returned too (least significant first).
operation RunModExp(cfg : Config, quantumExponent : Bool, exponentBits : Bool[]) : (Bool[], Bool[], Bool[]) {
    let size = Fold((acc, s) -> acc + s, 0, PieceSizes(cfg));
    let lens = PieceLowLens(cfg);
    use (qa, qb) = (Qubit[size], Qubit[size]);
    mutable acc = SplitRegister(cfg, qa);
    mutable work = SplitRegister(cfg, qb);
    InitRegister(cfg, acc, 1L);
    InitRegister(cfg, work, 0L);
    // The lookup output register is allocated after initialization, whose comparisons use a scratch
    // register as wide as the register being prepared.
    use qadd = Qubit[size];
    let addend = SplitRegister(cfg, qadd);

    let ne = Fold((a, s) -> a + s, 0, cfg.ExpWindows);
    mutable top = ne;
    mutable frequency = [];
    for wi in IndexRange(cfg.ExpWindows) {
        let c = cfg.ExpWindows[wi];
        set top -= c;
        use window = Qubit[c];
        if quantumExponent {
            ApplyToEach(H, window);
        } else {
            for b in 0..c - 1 {
                if exponentBits[top + b] {
                    X(window[b]);
                }
            }
        }
        let (ks, kinv) = cfg.Placeholder ? ([], []) | WindowConstants(cfg.Bases[wi], c, cfg.Modulus);
        // Every windowed multiplication of the same window size has the same structure, so when
        // counting it is traced once per window size and accounted for each time it occurs.
        if not cfg.Cache or BeginEstimateCaching($"windowed_multiplication_c{c}", 0) {
            MultiplyAdd(cfg, window, ks, acc, work, addend);
            MultiplyAdd(cfg, window, kinv, work, acc, addend);
            if cfg.Cache {
                EndEstimateCaching();
            }
        }
        set (acc, work) = (work, acc);
        if quantumExponent {
            set frequency = MeasureWindowSemiclassical(window, frequency, cfg.Placeholder);
        } else {
            ResetAll(window);
        }
    }
    let accBits = ResultArrayAsBoolArray(MResetEachZ(Flattened(acc)));
    let workBits = ResultArrayAsBoolArray(MResetEachZ(Flattened(work)));
    (accBits, workBits, frequency)
}

/// The same circuit, for counting with the QDK's logical resource counter.
operation CountModExp(cfg : Config) : Unit {
    let _ = RunModExp(cfg, true, []);
}

// ---------------------------------------------------------------------------------------------
// Pieces counted on their own (for the timing model and for checks)
// ---------------------------------------------------------------------------------------------

/// One lookup of `addressBits` address bits into n output qubits, without uncomputation.
operation CountLookup(numBits : Int, addressBits : Int) : Unit {
    use (address, out) = (Qubit[addressBits], Qubit[numBits]);
    Select([[false, size = numBits], size = 1 <<< addressBits], address, out);
    ResetAll(out);
}

/// One lookup and its measurement-based uncomputation.
operation CountLookupUnlookup(numBits : Int, addressBits : Int) : Unit {
    use (address, out) = (Qubit[addressBits], Qubit[numBits]);
    let data = [[false, size = numBits], size = 1 <<< addressBits];
    Select(data, address, out);
    Adjoint Select(data, address, out);
}

/// One Cuccaro addition on `len` qubits.
operation CountAdder(len : Int) : Unit {
    use (a, b) = (Qubit[len], Qubit[len]);
    CuccaroAdd(a, b);
}

/// Classical check of the adders: returns b + a (mod 2^n) and the carry-out variant's (sum, carry).
operation CheckAdders(n : Int, a : Int, b : Int) : (Int, Int, Bool) {
    use (qa, qb, qb2, z) = (Qubit[n], Qubit[n], Qubit[n], Qubit());
    ApplyXorInPlace(a, qa);
    ApplyXorInPlace(b, qb);
    ApplyXorInPlace(b, qb2);
    CuccaroAdd(qa, qb);
    CuccaroAddCarry(qa, qb2, z);
    let s1 = MeasureInteger(qb);
    let s2 = MeasureInteger(qb2);
    let carry = MResetZ(z) == One;
    // the addend must be restored
    Fact(MeasureInteger(qa) == a, "adder did not restore its addend");
    (s1, s2, carry)
}
