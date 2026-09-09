import os
import sys
import random

from sage.repl.load import load as load_sage_source

from utilities.sum_of_squares import sum_of_squares

# =============================================================================
# Load exactly one parameter set from the original, unmodified POKE_PKE.sage.
# =============================================================================

def load_poke_level(level=1, b_override=None):
    if level not in [1, 2, 3]:
        raise ValueError("level must be 1, 2, or 3")

    here = os.path.dirname(os.path.abspath(__file__))
    source_path = os.path.join(here, "POKE_PKE.sage")
    with open(source_path, "r") as f:
        source = f.read()

    marker = "    N = 100 # number of iterations"
    if marker not in source:
        raise ValueError("Could not locate the built-in benchmark marker")
    source = source.split(marker, 1)[0]

    lambdas = [128, 192, 256]
    old_loop = "for i, lambda_ in enumerate([128, 192, 256]):"
    new_loop = "for i, lambda_ in [(%d, %d)]:" % (level - 1, lambdas[level - 1])
    if old_loop not in source:
        raise ValueError("Could not locate the parameter loop")
    source = source.replace(old_loop, new_loop, 1)

    if b_override is not None:
        old_bs = "bs = [162, 243, 324]"
        if old_bs not in source:
            raise ValueError("Could not locate the POKE b parameter list")
        attack_bs = [162, 243, 324]
        attack_bs[level - 1] = int(b_override)
        source = source.replace(
            old_bs, "bs = %r" % attack_bs, 1
        )
        source = source.replace(
            "    assert B > 2**(2 * lambda_)",
            "    if B <= 2**(2 * lambda_):\n"
            "        print(\"WARNING: demonstration b does not meet the original security bound\")",
            1,
        )

    tmp = os.path.join(here, "._poke_level_%d.sage" % level)
    with open(tmp, "w") as f:
        f.write(source)
    try:
        load_sage_source(tmp, globals())
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


# =============================================================================
# Theta_dim4 external adapter.
# This imports the official Sage dimension-4 implementation without copying or
# modifying POKE_PKE.sage.  Set THETA_DIM4_PATH to the Theta_dim4_sage folder.
# =============================================================================

def load_theta4():
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.environ.get("THETA_DIM4_PATH"),
        os.path.join(here, "Theta_dim4_sage"),
        os.path.join(here, "Theta_dim4", "Theta_dim4_sage"),
        os.path.join(os.path.dirname(here), "Theta_dim4", "Theta_dim4_sage"),
    ]
    theta_path = None
    for path in candidates:
        if path and os.path.isdir(os.path.join(path, "pkg")):
            theta_path = os.path.abspath(path)
            break
    if theta_path is None:
        raise ImportError(
            "Theta_dim4_sage not found. Set THETA_DIM4_PATH to the directory "
            "containing pkg/, or run ./fetch_theta_dim4.sh"
        )
    if theta_path not in sys.path:
        sys.path.insert(0, theta_path)

    from pkg.isogenies.Kani_endomorphism import KaniEndoHalf
    from pkg.theta_structures.Tuple_point import TuplePoint
    from pkg.utilities.strategy import precompute_strategy_with_first_eval

    return {
        "path": theta_path,
        "KaniEndoHalf": KaniEndoHalf,
        "TuplePoint": TuplePoint,
        "strategy": precompute_strategy_with_first_eval,
    }


# =============================================================================
# Even-delta fault models.
# omega is odd.  We force omega_fault to remain odd and delta to be a nonzero
# even integer.  No change is made to any other POKE encryption operation.
# =============================================================================

def even_fault_scalar(omega, fault_type="even_random_byte", fault_bits=8,
                      bit_index=1):
    omega = ZZ(omega)

    if fault_type == "bit_flip":
        if int(bit_index) <= 0:
            raise ValueError("even delta requires bit_index >= 1")
        omega_fault = ZZ(int(omega) ^^ (1 << int(bit_index)))

    elif fault_type in ["even_random_byte", "even_random_low_bits"]:
        bits = 8 if fault_type == "even_random_byte" else int(fault_bits)
        if bits < 2:
            raise ValueError("even random replacement needs at least 2 bits")
        mask = (1 << bits) - 1
        old_low = int(omega) & mask
        # omega is odd.  Choosing a different odd replacement makes delta
        # nonzero and even, hence omega_fault remains invertible modulo 2^k.
        choices = [x for x in range(1, mask + 1, 2) if x != old_low]
        replacement = random.choice(choices)
        omega_fault = ZZ((int(omega) & ~mask) | replacement)

    else:
        raise ValueError("unknown even fault type")

    delta = ZZ(omega_fault - omega)
    assert delta != 0 and delta % 2 == 0
    assert omega_fault % 2 == 1
    return omega_fault, delta


def even_delta_candidates(fault_type="even_random_byte", fault_bits=8,
                          bit_index=1):
    if fault_type == "bit_flip":
        d = ZZ(1 << int(bit_index))
        return [-d, d]

    bits = 8 if fault_type == "even_random_byte" else int(fault_bits)
    bound = (1 << bits) - 1
    return [ZZ(d) for d in range(-bound, bound + 1)
            if d != 0 and d % 2 == 0]


# =============================================================================
# Faithful POKE encryption with random r3 and one external fault.
# The original POKE_PKE.sage is never edited.
# =============================================================================

def faithful_even_fault_encrypt(pkA, m, fault_type="even_random_byte",
                                fault_bits=8, bit_index=1):
    xP2, xQ2, xPQ2, xP3, xQ3, xPQ3, X_A, Y_A = pkA

    r3 = ZZ.random_element(B)             # Algorithm 3, Step 1
    d1, d2, d3, d4 = random_matrix(C)
    omega = random_unit(A)
    omega_inv = inverse_mod(omega, A)     # preserve supplied implementation
    omega_fault, delta = even_fault_scalar(
        omega, fault_type=fault_type,
        fault_bits=fault_bits, bit_index=bit_index
    )

    # psi : E0 -> EB, ker psi = <R0 + [r3]S0>
    _KB = _QB.ladder_3_pt(_PB, _PQB, r3)
    psi = KummerLineIsogeny(_E0, _KB, B)

    EB = psi.codomain()
    xP2_B = psi(xPA)
    xQ2_B = psi(xQA)
    X_B = psi(xX_0).curve_point()
    Y_B = psi(xY_0).curve_point()
    xXY_B = psi(xXY_0)
    if (X_B - Y_B)[0] != xXY_B.x():
        Y_B = -Y_B
    X_B, Y_B = d1*X_B + d2*Y_B, d3*X_B + d4*Y_B

    P2_B, Q2_B = lift_image_to_curve(PA, QA, xP2_B, xQ2_B, 4*A, B)
    P2_B = omega_fault * P2_B             # only injected fault
    Q2_B = omega_inv * Q2_B               # honest branch

    # psi' : EA -> EAB with the same r3
    EA = xP3.parent()
    xK = xQ3.ladder_3_pt(xP3, xPQ3, r3)
    psi_prime = KummerLineIsogeny(EA, xK, B)

    EAB = psi_prime.codomain().curve()
    xP2_AB = psi_prime(xP2)
    xQ2_AB = psi_prime(xQ2)
    xPQ2_AB = psi_prime(xPQ2)
    P2_AB = xP2_AB.curve_point()
    Q2_AB = xQ2_AB.curve_point()
    if (P2_AB - Q2_AB)[0] != xPQ2_AB.x():
        Q2_AB = -Q2_AB
    P2_AB *= omega
    Q2_AB *= omega_inv

    xX_AB = psi_prime(EA(X_A[0]))
    xY_AB = psi_prime(EA(Y_A[0]))
    xXY_AB = psi_prime(EA((X_A - Y_A)[0]))
    X_AB = xX_AB.curve_point()
    Y_AB = xY_AB.curve_point()
    if (X_AB - Y_AB)[0] != xXY_AB.x():
        Y_AB = -Y_AB
    X_AB, Y_AB = d1*X_AB + d2*Y_AB, d3*X_AB + d4*Y_AB

    ct_bytes = xof_encrypt(xof_kdf(X_AB[0], Y_AB[0]), m)
    ct = (EB.curve(), P2_B, Q2_B, X_B, Y_B,
          EAB, P2_AB, Q2_AB, ct_bytes)

    return ct


# =============================================================================
# Public scalar and torsion-image recovery.
# =============================================================================

def pairing_ratio(ct):
    EB, P_fault, Q_B, X_B, Y_B, EAB, P_AB, Q_AB, ct_bytes = ct

    # Project the implementation's 4A-torsion down to A-torsion, where the
    # masking relation using inverse_mod(omega,A) is exact.
    P_fault_A = 4 * P_fault
    Q_B_A = 4 * Q_B
    P0_A = 4 * PA
    Q0_A = 4 * QA

    lam_fault = P_fault_A.weil_pairing(Q_B_A, A)
    lam_base = P0_A.weil_pairing(Q0_A, A)
    exponent = discrete_log_pari(lam_fault, lam_base, A)

    return Zmod(A)(exponent) * Zmod(A)(B)^(-1)


def solve_linear_mod_power_of_two(delta, rhs, modulus):
    delta = ZZ(delta)
    rhs = ZZ(rhs)
    modulus = ZZ(modulus)
    g = gcd(delta, modulus)
    if rhs % g != 0:
        return []

    d_red = delta // g
    r_red = rhs // g
    m_red = modulus // g
    x0 = ZZ(Zmod(m_red)(r_red) * Zmod(m_red)(d_red)^(-1))
    return [ZZ(x0 + k*m_red) for k in range(ZZ(g))]


def recover_scalar_candidates(ct, candidates):
    ratio = pairing_ratio(ct)
    rhs = ZZ(ratio - 1)
    out = []
    seen = set()

    for delta in candidates:
        for omega_inv in solve_linear_mod_power_of_two(delta, rhs, A):
            if gcd(omega_inv, A) != 1:
                continue
            omega = ZZ(inverse_mod(omega_inv, A))
            omega_fault = ZZ(Zmod(A)(ratio) * Zmod(A)(omega))
            key = (ZZ(delta), omega, ZZ(omega_inv), omega_fault)
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "delta": ZZ(delta),
                "omega": omega,
                "omega_inv": ZZ(omega_inv),
                "omega_fault": omega_fault,
                "v2_delta": ZZ(delta).valuation(2),
            })
    return out


def recover_full_torsion_images(pkA, ct, candidate):
    xP2, xQ2, xPQ2, xP3, xQ3, xPQ3, X_A, Y_A = pkA
    EB, P_fault, Q_B, X_B, Y_B, EAB, P_AB, Q_AB, ct_bytes = ct
    modulus = 4*A

    omega = ZZ(candidate["omega"])
    omega_inv_A = ZZ(candidate["omega_inv"])
    omega_fault = ZZ(candidate["omega_fault"])

    if any(gcd(s, modulus) != 1 for s in [omega, omega_inv_A, omega_fault]):
        return None

    psiP = inverse_mod(omega_fault, modulus) * P_fault
    psiQ = inverse_mod(omega_inv_A, modulus) * Q_B
    psiP_prime = inverse_mod(omega, modulus) * P_AB
    psiQ_prime = inverse_mod(omega_inv_A, modulus) * Q_AB

    # Recover full source points on EA from the public x-only basis.
    P2 = xP2.curve_point()
    Q2 = xQ2.curve_point()
    if (P2 - Q2)[0] != xPQ2.x():
        Q2 = -Q2

    return {
        "P2": P2,
        "Q2": Q2,
        "psiP": psiP,
        "psiQ": psiQ,
        "psiP_prime": psiP_prime,
        "psiQ_prime": psiQ_prime,
    }


# =============================================================================
# Dimension-4 direct evaluator.  No HD inequality is checked here.
# =============================================================================

def find_kani_parameters(q, max_extra_e=64):
    """Find e,a1,a2 with q+a1^2+a2^2=2^e, without any torsion check."""
    q = ZZ(q)
    start = q.nbits()
    for e in range(start, start + int(max_extra_e) + 1):
        remainder = ZZ(2**e - q)
        if remainder <= 0:
            continue
        sol = sum_of_squares(remainder)
        if sol:
            a1, a2 = map(ZZ, sol)
            assert q + a1*a1 + a2*a2 == 2**e
            return e, a1, a2
    raise ValueError("No Kani parameters found in searched e-range")


def build_direct_embedding(theta4, P, Q, phiP, phiQ, q, e, a1, a2,
                           available_e):
    # Exactly the strategy setup used by Theta_dim4's SIDH attack, except that
    # we deliberately do not test whether available_e is theoretically enough.
    even_ai = a1 if a1 % 2 == 0 else a2
    m = ZZ(even_ai).valuation(2)
    f1 = ceil(ZZ(e) / 2)
    f2 = ZZ(e) - f1
    strategy1 = theta4["strategy"](f1, m, M=1, S=0.8, I=100)
    strategy2 = strategy1 if f2 == f1 else theta4["strategy"](
        f2, m, M=1, S=0.8, I=100
    )

    return theta4["KaniEndoHalf"](
        P, Q, phiP, phiQ, ZZ(q), ZZ(a1), ZZ(a2), ZZ(e),
        ZZ(available_e), strategy1, strategy2
    )


def eval_two_consistent(theta4, F, E1, E2, P, Q, component=2):
    TuplePoint = theta4["TuplePoint"]
    O1 = E1(0)
    O2 = E2(0)

    def evaluate(R):
        return F(TuplePoint(R, O1, O2, O2))[component]

    FP = evaluate(P)
    FQ = evaluate(Q)
    FPQ = evaluate(P - Q)

    # Resolve independent projective signs using the image of P-Q.
    choices = [
        (FP, FQ),
        (FP, -FQ),
        (-FP, FQ),
        (-FP, -FQ),
    ]
    for RP, RQ in choices:
        if RP - RQ == FPQ or RP - RQ == -FPQ:
            return RP, RQ
    raise ValueError("Could not make dimension-4 evaluations sign-consistent")
