"""Focused Level-5 POKE fault attack with reverse-5 recovery of r3."""

import os
import sys
import time
from sage.repl.load import load as load_sage_source

from montgomery_isogenies.kummer_line import KummerLine
from montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
from utilities.supersingular import torsion_basis


DEPTH = 5


def show_progress(done, total, started):
    """Draw an in-place progress bar with elapsed time and worst-case ETA."""
    elapsed = time.time() - started
    ratio = float(done) / total if total else 1.0
    width = 30
    filled = min(width, int(width * ratio))
    bar = "#" * filled + "-" * (width - filled)
    eta = elapsed * (total - done) / done if done else 0.0
    print("\r    [%s] %6.2f%%  %d/%d  elapsed %5.1f min  ETA(max) %5.1f min" %
          (bar, 100.0 * ratio, done, total, elapsed / 60.0, eta / 60.0),
          end="", flush=True)


def load_helpers():
    scope = dict(globals())
    scope["__name__"] = "faulty_poke_helpers"
    load_sage_source(os.path.join(os.path.dirname(__file__),
                                  "faulty_poke_direct.sage"), scope)
    return scope


def reverse_chains(EB):
    """Enumerate the 4*3^4 non-backtracking paths leaving EB."""
    states = [(KummerLine(EB), None, ())]
    for _ in range(DEPTH):
        children = []
        for line, previous_j, chain in states:
            curve = line.curve()
            P, Q = torsion_basis(curve, 3)
            for K in (P, Q, P + Q, P - Q):
                phi = KummerLineIsogeny(line, line(K[0]), 3)
                next_line = phi.codomain()
                if (previous_j is not None and
                    next_line.curve().j_invariant() == previous_j):
                    continue
                children.append((next_line, curve.j_invariant(),
                                 chain + (phi,)))
        states = children
    return [(line, chain) for line, _, chain in states]


def push_x(chain, point):
    for phi in chain:
        point = phi(point)
    return point


def prefix_images(H, images, chain, A):
    P, Q = images["psiP"], images["psiQ"]
    _, xP, xQ, xPQ = H["point_to_xonly"](P, Q)
    xP, xQ, xPQ = (push_x(chain, T) for T in (xP, xQ, xPQ))
    P, Q = xP.curve_point(), xQ.curve_point()
    if (P - Q)[0] != xPQ.x():
        Q = -Q
    inv = ZZ(inverse_mod(3**DEPTH, 4*A))
    return inv*P, inv*Q


def recover_r(H, theta4, F, Epre, EB, q):
    """Recover r mod q by DLP, then its five remaining base-3 digits."""
    from pkg.utilities.discrete_log import ell_discrete_log_pari

    P, Q = 3**DEPTH*H["PB"], 3**DEPTH*H["QB"]
    U, V = H["eval_two_consistent"](
        theta4, F, H["E0"], Epre, P, Q, component=2
    )
    r0 = ZZ(ell_discrete_log_pari(Epre, -U, V, q)) % q

    U, V = H["eval_two_consistent"](
        theta4, F, H["E0"], Epre, H["PB"], H["QB"], component=2
    )
    line = KummerLine(Epre)
    answers = []
    for high in range(3**DEPTH):
        r = ZZ(r0 + high*q)
        K = U + r*V
        if K.is_zero() or 3**(DEPTH - 1)*K == Epre(0):
            continue
        tail = KummerLineIsogeny(line, line(K[0]), 3**DEPTH)
        if tail.codomain().curve().j_invariant() == EB.j_invariant():
            answers.append(r)
    return answers


def decrypt_with_r(H, pk, ct, r):
    """Rebuild psi and psi' from r, recover D, and decrypt."""
    xP2, xQ2, xPQ2, xP3, xQ3, xPQ3, X_A, Y_A = pk
    EB, _, _, X_B, Y_B, EAB, _, _, ct_bytes = ct

    K = H["_QB"].ladder_3_pt(H["_PB"], H["_PQB"], r)
    psi = KummerLineIsogeny(H["_E0"], K, H["B"])
    if psi.codomain().curve().j_invariant() != EB.j_invariant():
        raise ValueError("recovered r3 does not reproduce public EB")
    X, Y = psi(H["xX_0"]).curve_point(), psi(H["xY_0"]).curve_point()
    if (X - Y)[0] != psi(H["xXY_0"]).x():
        Y = -Y
    wp = X.weil_pairing(Y, H["C"])
    d1, d2 = H["BiDLP"](X_B, X, Y, H["C"], ePQ=wp)
    d3, d4 = H["BiDLP"](Y_B, X, Y, H["C"], ePQ=wp)

    EA = xP3.parent()
    K = xQ3.ladder_3_pt(xP3, xPQ3, r)
    psi_prime = KummerLineIsogeny(EA, K, H["B"])
    if psi_prime.codomain().curve().j_invariant() != EAB.j_invariant():
        raise ValueError("recovered r3 does not reproduce public EAB")
    X = psi_prime(EA(X_A[0])).curve_point()
    Y = psi_prime(EA(Y_A[0])).curve_point()
    if (X - Y)[0] != psi_prime(EA((X_A - Y_A)[0])).x():
        Y = -Y
    X, Y = d1*X + d2*Y, d3*X + d4*Y
    return H["xof_encrypt"](H["xof_kdf"](X[0], Y[0]), ct_bytes)


def attack(seed=1):
    t0 = time.time()
    print("[1/6] Loading original Level-5 parameters and Theta_dim4...",
          flush=True)
    set_random_seed(seed)
    H = load_helpers()
    H["set_random_seed"](seed)
    H["random"].seed(seed)
    H["load_poke_level"](3)               # original b=324
    theta4 = H["load_theta4"]()

    print("[2/6] Generating a key and one faulted ciphertext...", flush=True)
    _, pk = H["keygenA"]()
    ct = H["faithful_even_fault_encrypt"](
        pk, os.urandom(32), fault_type="bit_flip", bit_index=1
    )
    print("[3/6] Recovering public masking-scalar candidates...", flush=True)
    scalars = H["recover_scalar_candidates"](
        ct, H["even_delta_candidates"]("bit_flip", 8, 1)
    )
    print("      recovered %d masking candidates" % len(scalars), flush=True)
    print("[4/6] Enumerating 324 reverse paths and preparing Kani parameters...",
          flush=True)
    chains = reverse_chains(ct[0])

    q = ZZ(3**(H["b"] - DEPTH))
    e, a1, a2 = H["find_kani_parameters"](q)
    attempts = 0
    checked = 0
    recovered_r = None
    total_candidates = len(scalars) * len(chains)
    search_started = time.time()
    print("[5/6] Testing at most %d public HD candidates..." % total_candidates,
          flush=True)

    for scalar in scalars:
        images = H["recover_full_torsion_images"](pk, ct, scalar)
        if images is None:
            checked += len(chains)
            show_progress(checked, total_candidates, search_started)
            continue
        for Epre_line, chain in chains:
            attempts += 1
            checked += 1
            answers = []
            try:
                P, Q = prefix_images(H, images, chain, H["A"])
                F = H["build_direct_embedding"](
                    theta4, H["PA"], H["QA"], P, Q, q, e, a1, a2,
                    H["a"] + 2
                )
                # Public validity filter: wrong paths fail during gluing.
                H["eval_two_consistent"](
                    theta4, F, H["E0"], Epre_line.curve(),
                    H["X0"], H["Y0"], component=2
                )
                answers = recover_r(H, theta4, F, Epre_line.curve(), ct[0], q)
            except (AssertionError, ValueError, ZeroDivisionError):
                pass
            show_progress(checked, total_candidates, search_started)
            if len(answers) == 1:
                recovered_r = answers[0]
                break
        if recovered_r is not None:
            break

    print(flush=True)
    if recovered_r is None:
        raise ValueError("No public candidate recovered r3")
    print("      unique public candidate found after %d HD attempts" % attempts,
          flush=True)
    print("[6/6] Validating public EB/EAB and decrypting...", flush=True)
    plaintext = decrypt_with_r(H, pk, ct, recovered_r)
    elapsed = time.time() - t0

    print("masking candidates:", len(scalars))
    print("reverse paths:", len(chains))
    print("HD candidates attempted:", attempts)
    print("public EB/EAB validation: passed")
    print("recovered r3:", recovered_r)
    print("recovered plaintext (hex):", plaintext.hex())
    print("total time: %.3f s" % elapsed)
    return recovered_r, plaintext


if __name__ == "__main__":
    attack(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
