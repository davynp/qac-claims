<!-- Rendered by `qacclaims build` from narrative.md. Edit narrative.md; this file is regenerated. -->

# Reproducing the RSA-2048 resource estimate of Gidney and Ekerå (2021)

*Siya Mehta, PlexTech, UC Berkeley.*
*Numbers in this page are rendered from the claim record; the record's numbers come from the executed notebook in `reproductions/qsharp-modexp-qre-2026-11/`.*

## The claim

[Gidney and Ekerå (2021)](https://doi.org/10.22331/q-2021-04-15-433) estimate what it would take to break
RSA-2048, the most widely used RSA key size, with Shor-type
algorithms on an error-corrected quantum computer. Read literally, the paper claims:

> A 2048-bit RSA integer could be factored in about 8 hours with about 20 million noisy qubits.

Two layers of numbers sit behind that sentence.

**Logical counts.** Almost all of the work is one modular exponentiation. The paper makes it cheap by
combining several techniques: Ekerå and Håstad's variant of Shor's algorithm, which needs fewer
multiplications; windowed arithmetic, which turns many controlled additions into a few table lookups
followed by additions; the coset representation of modular integers, which lets cheap non-modular adders
do modular arithmetic approximately; and oblivious carry runways, which split each register into pieces
that are added in parallel. The paper reports about 2.7 × 10⁹
Toffoli gates per run (Table 1). For the number of logical qubits it gives only an expression,
$3n + 0.002\,n \lg n$, and no value at this size, so this record does not compare logical qubits.

**Physical resources.** The paper assumes a planar grid of qubits with nearest-neighbour connectivity, a
physical gate error rate of 1 × 10⁻³, a surface-code
cycle of 1 × 10⁻⁶ s, and a reaction time of
1 × 10⁻⁵ s for the control system to decode a measurement and
act on it. It lays the computation out with Gidney–Fowler CCZ factories feeding ripple-carry additions
that are limited by that reaction time, uses distance-27
surface-code patches, and accepts a retry risk of 0.31 per run
(Table 3). Including repeated runs, it arrives at about 20 million
physical qubits and about 8 hours (title).

## What QAC did

We built the modular exponentiation of the paper's Section 2 as a circuit in Q#, the programming language of
the [QDK](https://github.com/microsoft/qdk) (version 1.33.1), and
costed it ourselves rather than evaluating the paper's formulas. The
[executed notebook](https://github.com/QuantumApplicationsCenter/qac-claims/blob/main/claims/gidney-ekera-rsa2021/reproductions/qsharp-modexp-qre-2026-11/notebook.ipynb)
walks through every step.

**The circuit.** Each register is stored in runway pieces with coset padding at the top. Each windowed
multiplication is two multiply-adds and a swap. Each multiply-add looks up a precomputed multiple of a
constant, addressed by a window of exponent bits and a window of the other register's bits, adds it into
every piece of the target in parallel with Cuccaro's ripple-carry adder, and uncomputes the lookup by
measurement. The exponent is read out by the semi-classical Fourier transform, a few qubits at a time.

**Checking it.** On small instances the notebook checks the adders on random inputs, and the whole
exponentiation exactly against modular powers, with one, two and three runway pieces. It checks that
preparing the padding and runways in superposition fails exactly as often as an exhaustive classical count
predicts; with tiny paddings that is often, which is the approximation the paper bounds. And it runs order
finding end to end: the measured frequencies land on the peaks of the ideal distribution, so the arithmetic
keeps the phases that Shor's algorithm needs. Building it exposed two subtle mistakes of ours, an inverted
comparison in the initialization and an ordering of the runway setup that made overflow more likely; the
checks caught both, and both are fixed.

**Logical counts.** With Table 3's parameters ($n_e$ = 3,029
exponent qubits, windows of 5 and
5 bits, runways every
1,024 bits, padding of
38 bits), the circuit needs
2.69 × 10⁹ Toffoli gates per run, a ratio of
1.00× to the paper's figure, and
6,430 logical qubits at its peak. As a cross-check, the
notebook also downloads the authors' own cost script from arXiv and reruns it: its parameter search selects
exactly Table 3's configuration, and its Toffoli count agrees with the circuit's to within a fraction of a
percent.

**Physical resources.** We rebuilt the paper's layout and timing model (Secs. 2.10–2.14) as code, with each
constant tied to its place in the text, and applied it to the circuit's own counts: lookups limited by code
depth, additions limited by the reaction time, runway pieces in parallel, and failures from the surface code,
the factories, the approximation and post-processing combined into a retry risk. This gives
22.33 million physical qubits
(22,325,184) and
5.1 hours per run, or
6.8 hours expected including repeated runs. An independent
estimate with the QDK resource estimator, using magic state cultivation and an 8T-to-CCX factory in place of
the paper's CCZ factories, chooses distance-27 patches on
its own and needs 19.1 million physical qubits, but
2.7 days of expected runtime.

## How to read the comparison

The Toffoli count compares the same construction at the same parameters, and it agrees closely
(1.00×). The physical figures use the paper's own model, so their
ratios (1.12× for physical qubits and
0.85× for runtime) mostly reflect what the text does and does
not specify. The text gives the factory footprint only for a larger level-1 distance than Table 3 uses,
and its own hand calculation at that layout (Sec. 2.14) also lands above the headline; the paper says the
headline comes from its script, which finds a more compact layout at Table 3's level-1 distance of
15, a layout the text does not describe. The runtime is a little shorter than the rounded headline.

The estimator's figures show which assumptions carry the claim. Its qubit count and code distance land
close to the paper's, but its runtime is far longer, because it executes the Toffolis one after another and
has no notion of the reaction-limited additions that run in parallel runway pieces. The paper's
8 hours is a property of that architecture, not just of the gate
count.

## What this does not show, and what comes next

- It does not reproduce the title's qubit count from the text alone: the more compact factory layout behind
  it is in the authors' script, not in the paper.
- The distillation error per CCZ state is taken from the paper's aggregate figure (Sec. 2.12), not from a
  simulation of the factory, and the CCZ factory itself is not rebuilt.
- The full-size circuit is counted, not simulated; the order-finding check uses a small modulus.
- No logical-qubit comparison is made, because the paper prints only an asymptotic expression.
- A natural next step is to rebuild the factory layout of the paper's Refs. 36 and 37 at the level-1
  distance Table 3 uses, which would make the physical-qubit comparison exact.

## Reproducing this result

The [notebook](https://github.com/QuantumApplicationsCenter/qac-claims/blob/main/claims/gidney-ekera-rsa2021/reproductions/qsharp-modexp-qre-2026-11/notebook.ipynb)
runs in about 209 seconds on a laptop, in the environment pinned by
`requirements.lock.txt` (Python 3.10.18). It downloads the authors'
script from arXiv, so it needs network access. From a checkout of this repository:

```bash
python3.10 -m venv .venv && source .venv/bin/activate
pip install -r claims/gidney-ekera-rsa2021/reproductions/qsharp-modexp-qre-2026-11/requirements.lock.txt
pip install -e ".[execute]"
qacclaims execute claims/gidney-ekera-rsa2021/reproductions/qsharp-modexp-qre-2026-11
```

The last command re-executes the notebook and reports any difference from the committed `results.json`.
