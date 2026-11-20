"""The physical layout and timing model of Gidney and Ekerå (2021), Section 2.10-2.14, as code.

Every constant below is quoted from the paper's text, with its locator. Page numbers follow the
published version (Quantum 5, 433), which has the same pagination as arXiv:1905.09749v3. The model is
applied to counts taken from the Q# circuit, not to the paper's formulas.

The text gives the CCZ factory footprint only for a level-1 code distance of 17 (Sec. 2.10, p. 11),
although Table 3's n = 2048 row uses 15. The rebuilt layout therefore uses 17, which is the layout
Sec. 2.14 (p. 13) evaluates; the paper says there that its published estimate comes from a script that
can use 15, which "adjusts the layout in several fortuitous ways".
"""

from __future__ import annotations

import math

# --- Physical assumptions (abstract, p. 1; Table 2, p. 3) -------------------------------------
PHYSICAL_ERROR_RATE = 1e-3
CYCLE_TIME_S = 1e-6
REACTION_TIME_S = 10e-6

# --- Code distances ------------------------------------------------------------------------------
CODE_DISTANCE = 27            # level-2 distance, Table 3 (n = 2048 row, p. 16); Sec. 2.10-2.14
LEVEL1_DISTANCE = 17          # the level-1 distance for which Sec. 2.10 (p. 11) states the footprint


def physical_qubits_per_patch(d: int) -> int:
    """2(d + 1)^2 physical qubits per logical qubit (Sec. 1.1, p. 2; Sec. 2.14, p. 13; Fig. 8, p. 14)."""
    return 2 * (d + 1) ** 2


def logical_error_per_patch_cycle(d: int) -> float:
    """10^-ceil(d/2 + 1) per logical qubit per cycle at a 10^-3 gate error rate (Sec. 2.13, p. 13)."""
    return 10.0 ** (-math.ceil(d / 2 + 1))


# --- Layout of one runway piece (Sec. 2.10, p. 11) ---------------------------------------------
FACTORIES_PER_PIECE = 14      # CCZ factories needed to run at the reaction-limited rate at d = 27
FACTORY_WIDTH, FACTORY_HEIGHT = 15, 8   # footprint at level-1 distance 17
FACTORY_ROWS = 2              # "laid out into 2 rows of 7, with single logical qubit gaps"
OPERATING_ROWS = (
    FACTORY_ROWS * FACTORY_HEIGHT  # "2 * 8 for the CCZ factories"
    + 3                            # "three for the ripple-carry operating area"
    + 6                            # "six for AutoCCZ fixup boxes"
    + 8                            # "eight for routing data qubits"
)                                  # = 33
REGISTERS_IN_LAYOUT = 3        # target, input and lookup-output registers (Sec. 2.10, Figs. 4-6)


PIECE_WIDTH = 113             # "The total width of a piece is 15 * 7 + 7 = 113 logical qubits" (Sec. 2.10, p. 11)


def piece_width() -> int:
    """The stated 113 logical qubits. The paper's own expression, 15 * 7 + 7, is 112; the stated
    value is used (it matches one more gap column, as in the paper's script)."""
    per_row = FACTORIES_PER_PIECE // FACTORY_ROWS
    assert FACTORY_WIDTH * per_row + per_row == PIECE_WIDTH - 1
    return PIECE_WIDTH


def layout(*, pieces: int, piece_qubits: int, d: int = CODE_DISTANCE) -> dict:
    """Board size for `pieces` runway pieces whose registers hold `piece_qubits` qubits per piece.

    Each register's share of a piece (c_sep + c_pad qubits) is packed into rows of the piece's width;
    the paper's "another 30 rows (approximately)" is 3 registers of ceil(1062 / 113) = 10 rows.
    """
    width = piece_width()
    register_rows = REGISTERS_IN_LAYOUT * math.ceil(piece_qubits / width)
    height = OPERATING_ROWS + register_rows
    patches = width * pieces * height
    factory_patches = FACTORIES_PER_PIECE * FACTORY_WIDTH * FACTORY_HEIGHT * pieces
    per_patch = physical_qubits_per_patch(d)
    return {
        "board_width": width * pieces,
        "board_height": height,
        "logical_patches": patches,
        "factory_patches": factory_patches,
        "physical_qubits": patches * per_patch,
        "factory_qubits": factory_patches * per_patch,
        "compute_qubits": (patches - factory_patches) * per_patch,
        "factories": FACTORIES_PER_PIECE * pieces,
    }


# --- Timing of one lookup addition (Sec. 2.11, p. 11) -------------------------------------------
OTHER_LOOKUP_ADDITION_TIME_S = 1e-3   # "uncomputing the lookup and rearranging the rows ... 1 millisecond"


def lookup_time_s(lookup_toffolis: int, d: int = CODE_DISTANCE) -> float:
    """Code-depth limited lookup with double-access hallways: cycle * d/2 per sequential Toffoli
    (the paper writes 1 us * d/2 * 2^(c_exp + c_mul))."""
    return CYCLE_TIME_S * d / 2 * lookup_toffolis


def addition_time_s(adder_toffolis_per_piece: int) -> float:
    """Reaction-limited ripple-carry addition: one reaction time per sequential Toffoli of a piece's
    adder (the paper writes 2(c_sep + c_pad) * 10 us). Pieces run in parallel."""
    return REACTION_TIME_S * adder_toffolis_per_piece


def runway_reduction_time_s(reduction_toffolis: int) -> float:
    """The Fig. 3 runway reductions are ripple-carry additions too; reaction limited."""
    return REACTION_TIME_S * reduction_toffolis


# --- Errors (Secs. 2.9, 2.12, 2.13, 2.1, 3.1) ---------------------------------------------------
# Sec. 2.12 (p. 13): "approximately 3 billion CCZ states ... a total distillation error of 6.4%"
# at level-1 distance 17 and level-2 distance 27. The text gives no per-state figure; the model uses
# the ratio, scaled to the circuit's own Toffoli count.
DISTILLATION_ERROR_REPORTED = 0.064
DISTILLATION_STATES_REPORTED = 3e9
POST_PROCESSING_FAILURE = 0.01   # "probability above 99% of recovering d" (Sec. 2.1, p. 5; Sec. 3.2)


def topological_error(non_factory_patches: int, cycles: float, d: int = CODE_DISTANCE) -> float:
    """Sec. 2.13 (p. 13): per-patch-cycle error times the non-distillation patches times the cycles."""
    return logical_error_per_patch_cycle(d) * non_factory_patches * cycles


def distillation_error(toffolis: int) -> float:
    return toffolis * DISTILLATION_ERROR_REPORTED / DISTILLATION_STATES_REPORTED


def approximation_error(lookup_additions: int, n: int, runway_sep: int, padding: int) -> float:
    """Sec. 2.9 (pp. 9-10): deviation per addition at most n / (c_sep 2^c_pad), subadditive; the trace
    distance is at most 2 sqrt(deviation)."""
    deviation = lookup_additions * n / (runway_sep * 2 ** padding)
    return 2 * math.sqrt(deviation)


def union(*ps: float) -> float:
    """Combine error probabilities assuming independence (Sec. 3.1, p. 15)."""
    keep = 1.0
    for p in ps:
        keep *= 1 - min(p, 1.0)
    return 1 - keep
