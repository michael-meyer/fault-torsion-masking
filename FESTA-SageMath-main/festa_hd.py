"""Dimension-4 Kani/HD adapter used by the FESTA fault-attack PoC.

The code mirrors the direct evaluator used by the POKE Level-5 attack in this
artifact.  It consumes only an isogeny degree and its action on a sufficiently
large 2-power torsion basis.  In particular it does not need the FESTA secret
key or the kernel scalar of the isogeny being represented.
"""

import os
import sys

from sage.all import ZZ, ceil

from utilities.sum_of_squares import sum_of_squares


def load_theta4():
    """Load the bundled Theta_dim4 implementation."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.environ.get("THETA_DIM4_PATH"),
        os.path.join(here, "Theta_dim4_sage"),
        # Reuse exactly the same bundled implementation as the POKE attack.
        os.path.join(
            os.path.dirname(here), "POKE-Level5-Fault-Attack", "Theta_dim4_sage"
        ),
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
            "containing pkg/."
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


def find_kani_parameters(q, max_extra_e=96):
    """Find e,a1,a2 such that q+a1^2+a2^2 = 2^e."""
    q = ZZ(q)
    start = q.nbits()
    for e in range(start, start + int(max_extra_e) + 1):
        remainder = ZZ(2**e - q)
        if remainder <= 0:
            continue
        sol = sum_of_squares(remainder)
        if sol:
            a1, a2 = map(ZZ, sol)
            assert q + a1 * a1 + a2 * a2 == 2**e
            return ZZ(e), a1, a2
    raise ValueError("No Kani parameters found in searched e-range")


def build_direct_embedding(theta4, P, Q, phiP, phiQ, degree, available_e,
                           kani_parameters=None):
    """Build a KaniEndoHalf representation of an odd-degree isogeny.

    P,Q are a basis of E[2^available_e], and phiP,phiQ are their images under
    the unknown ``degree``-isogeny.  This is exactly the information recovered
    from the faulted FESTA masking.
    """
    if kani_parameters is None:
        e, a1, a2 = find_kani_parameters(degree)
    else:
        e, a1, a2 = kani_parameters
        e, a1, a2 = ZZ(e), ZZ(a1), ZZ(a2)

    if ZZ(available_e) < ceil(ZZ(e) / 2) + 2:
        raise ValueError(
            "insufficient 2-power torsion for KaniEndoHalf: have 2^%s, need "
            "at least 2^%s" % (available_e, ceil(ZZ(e) / 2) + 2)
        )

    even_ai = a1 if a1 % 2 == 0 else a2
    m = ZZ(even_ai).valuation(2)
    f1 = ceil(ZZ(e) / 2)
    f2 = ZZ(e) - f1
    strategy1 = theta4["strategy"](f1, m, M=1, S=0.8, I=100)
    strategy2 = (
        strategy1
        if f2 == f1
        else theta4["strategy"](f2, m, M=1, S=0.8, I=100)
    )

    F = theta4["KaniEndoHalf"](
        P,
        Q,
        phiP,
        phiQ,
        ZZ(degree),
        a1,
        a2,
        e,
        ZZ(available_e),
        strategy1,
        strategy2,
    )
    return F, (e, a1, a2)


def _make_pair_sign_consistent(FP, FQ, FPQ):
    """Resolve the independent projective signs using the image of P-Q."""
    for RP, RQ in ((FP, FQ), (FP, -FQ), (-FP, FQ), (-FP, -FQ)):
        if RP - RQ == FPQ or RP - RQ == -FPQ:
            return RP, RQ
    raise ValueError("Could not make dimension-4 evaluations sign-consistent")


def eval_forward_two_consistent(theta4, F, E_source, E_target, P, Q):
    """Evaluate the represented isogeny on P,Q, up to one overall sign."""
    TuplePoint = theta4["TuplePoint"]
    O1 = E_source(0)
    O2 = E_target(0)

    def evaluate(R):
        # For F=[[alpha_1, hat(sigma)],[-sigma, alpha_2]], component 2 is
        # -sigma(R) on input (R,0,0,0).  The common minus sign is harmless.
        return F(TuplePoint(R, O1, O2, O2))[2]

    FP = evaluate(P)
    FQ = evaluate(Q)
    FPQ = evaluate(P - Q)
    return _make_pair_sign_consistent(FP, FQ, FPQ)


def eval_dual_two_consistent(theta4, F, E_source, E_target, P, Q):
    """Evaluate the dual isogeny on target points P,Q.

    In the same Kani block matrix, input (0,0,R,0) has first component
    hat(sigma)(R).  This is the operation needed by FESTA's
    ComputeCanonicalKernel step.
    """
    TuplePoint = theta4["TuplePoint"]
    O1 = E_source(0)
    O2 = E_target(0)

    def evaluate(R):
        return F(TuplePoint(O1, O1, R, O2))[0]

    FP = evaluate(P)
    FQ = evaluate(Q)
    FPQ = evaluate(P - Q)
    return _make_pair_sign_consistent(FP, FQ, FPQ)
