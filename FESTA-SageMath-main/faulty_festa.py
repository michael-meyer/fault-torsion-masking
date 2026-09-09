"""FESTA message-recovery fault attack following Algorithm 7."""

import random
import sys
import time

from sage.all import Matrix, ZZ, Zmod, gcd, inverse_mod, set_random_seed

from festa import FESTA
from festa_hd import (
    build_direct_embedding,
    eval_dual_two_consistent,
    find_kani_parameters,
    load_theta4,
)
from parameters.params import (
    festa_params_128,
    festa_params_192,
    festa_params_256,
    festa_params_toy,
)
from utilities.discrete_log import windowed_pohlig_hellman
from utilities.masking import canonical_matrix
from utilities.pairing import weil_pairing_pari
from utilities.supersingular import compute_canonical_kernel, torsion_basis
from montgomery_isogenies.kummer_line import KummerLine
from montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
from montgomery_isogenies.isogenies_x_only import evaluate_isogeny_x_only


PARAMETERS = {
    "toy": festa_params_toy,
    "128": festa_params_128,
    "192": festa_params_192,
    "256": festa_params_256,
}


# =============================================================================
# Fault model.
# For simplicity, to show the recoverability, we let Delta be even here.
# =============================================================================

def bit_flip_delta_candidates(bit_index=1):
    """Return the Delta candidates for one bit-flip fault."""
    bit_index = int(bit_index)
    if bit_index <= 0:
        raise ValueError("use bit_index >= 1 so the faulted scalar stays odd")
    d = ZZ(1 << bit_index)
    return (-d, d)


def solve_linear_mod_power_of_two(delta, rhs, modulus):
    """Solve Delta*x = rhs modulo 2^b."""
    delta = ZZ(delta)
    rhs = ZZ(rhs)
    modulus = ZZ(modulus)
    g = gcd(delta, modulus)
    if rhs % g != 0:
        return []

    d_red = delta // g
    r_red = rhs // g
    m_red = modulus // g
    x0 = ZZ(Zmod(m_red)(r_red) * Zmod(m_red)(d_red) ** (-1))
    return [ZZ(x0 + k * m_red) for k in range(int(g))]


def _simultaneous_sign_key(alpha, beta, alpha_fault, modulus):
    """Canonical key modulo <-I_2>."""
    modulus = ZZ(modulus)
    raw = tuple(ZZ(x) % modulus for x in (alpha, beta, alpha_fault))
    neg = tuple((-x) % modulus for x in raw)
    return min(raw, neg)


def _matrix_key(M):
    M = canonical_matrix(M)
    return tuple(ZZ(x) for x in M.list())


def recover_secret_scalars(ctx, pk_compressed, ct, deltas):
    """Recover alpha', beta, and alpha as in Algorithm 7, Steps 1--4."""
    EA, PA, QA = ctx.decompress_public_key(pk_compressed)
    E1, R1_fault, S1, E2, R2, S2 = ct

    N = ZZ(ctx.l_power)
    R = Zmod(N)

    # e(R1',S1) = e(P0,Q0)^(alpha' beta d1)
    base1 = weil_pairing_pari(ctx.Pb, ctx.Qb, N)
    pair1 = weil_pairing_pari(R1_fault, S1, N)
    alpha_fault_beta_d1 = windowed_pohlig_hellman(
        pair1, base1, int(ctx.b), ctx.window
    )

    # e(R2,S2) = e(PA,QA)^(alpha beta d2)
    base2 = weil_pairing_pari(PA, QA, N)
    pair2 = weil_pairing_pari(R2, S2, N)
    alpha_beta_d2 = windowed_pohlig_hellman(
        pair2, base2, int(ctx.b), ctx.window
    )

    alpha_fault_beta = R(alpha_fault_beta_d1) * R(ctx.d1) ** (-1)
    alpha_beta = R(alpha_beta_d2) * R(ctx.d2) ** (-1)
    rhs = ZZ(alpha_fault_beta - alpha_beta)

    out = []
    seen = set()
    for delta in deltas:
        for beta_int in solve_linear_mod_power_of_two(delta, rhs, N):
            beta = R(beta_int)
            if gcd(ZZ(beta), N) != 1:
                continue

            alpha = alpha_beta * beta ** (-1)
            alpha_fault = alpha + R(delta)
            if gcd(ZZ(alpha), N) != 1 or gcd(ZZ(alpha_fault), N) != 1:
                continue

            key = _simultaneous_sign_key(alpha, beta, alpha_fault, N)
            if key in seen:
                continue
            seen.add(key)

            T_raw = Matrix(R, 2, 2, [alpha, 0, 0, beta])
            out.append(
                {
                    "delta": ZZ(delta),
                    "alpha": alpha,
                    "beta": beta,
                    "alpha_fault": alpha_fault,
                    "T": canonical_matrix(T_raw),
                    "sign_key": key,
                }
            )

    return (EA, PA, QA), out


def recover_torsion_images(candidate, ct):
    """Recover the unmasked torsion points as in Algorithm 7, Step 5."""
    E1, R1_fault, S1, E2, R2, S2 = ct
    alpha = candidate["alpha"]
    beta = candidate["beta"]
    alpha_fault = candidate["alpha_fault"]

    return {
        "psi1_P": alpha_fault ** (-1) * R1_fault,
        "psi1_Q": beta ** (-1) * S1,
        "psi2_P": alpha ** (-1) * R2,
        "psi2_Q": beta ** (-1) * S2,
    }


def recover_message_from_preimage(ctx, s, t, T):
    """Complete message recovery as in Algorithm 7, Steps 10--11."""
    s = Zmod(ctx.d1)(s)
    t = Zmod(ctx.d2)(t)
    x, X = ctx.H(s)
    r = Zmod(ctx.d2)(t - x)
    R = canonical_matrix(X.inverse() * T)
    m_prime = Zmod(ctx.d1)(s - ctx.G(r, R))
    m_prime_int = ZZ(m_prime)
    padding_ok = (m_prime_int & ((1 << int(ctx.k)) - 1)) == 0
    m = m_prime_int >> int(ctx.k)
    return {
        "s": ZZ(s),
        "t": ZZ(t),
        "T": canonical_matrix(T),
        "r": ZZ(r),
        "R": R,
        "m_prime": m_prime_int,
        "padding_ok": padding_ok,
        "message": ZZ(m),
    }


def prepare_psi1_d1_over_3_paths(ctx, E1):
    """Enumerate the four reverse 3-isogenies used for the d1/3 HD representation."""
    if ZZ(ctx.d1) % 3 != 0:
        raise ValueError("psi1 degree-d1/3 representation expects 3 | d1")

    q1 = ZZ(ctx.d1) // 3
    P3, Q3 = torsion_basis(E1, ZZ(3), even_power=ctx.b)
    kernels = (P3, Q3, P3 + Q3, P3 - Q3)
    L1 = KummerLine(E1)
    p = E1.base_ring().characteristic()

    paths = []
    for i, K in enumerate(kernels):
        rho = KummerLineIsogeny(L1, L1(K), ZZ(3))
        Epre = rho.codomain().curve()
        Epre.set_order((p + 1) ** 2, num_checks=0)
        target_basis = torsion_basis(Epre, q1, even_power=ctx.b)
        paths.append(
            {
                "path_index": i,
                "rho": rho,
                "Epre": Epre,
                "target_basis": target_basis,
            }
        )
    return q1, paths


def psi1_d1_over_3_images(ctx, images, path):
    """Recover the degree-d1/3 images from the candidate psi1 images."""
    N = ZZ(ctx.l_power)
    rhoP, rhoQ = evaluate_isogeny_x_only(
        path["rho"], images["psi1_P"], images["psi1_Q"], N, ZZ(3)
    )
    inv3 = ZZ(inverse_mod(3, N))
    return inv3 * rhoP, inv3 * rhoQ


def complete_public_candidate(
    ctx,
    theta4,
    public_key,
    ct,
    candidate,
    kani1_d1_over_3,
    kani2,
    psi1_degree_d1_over_3,
    psi1_paths,
    source_basis_d1_over_3,
    source_basis_d2,
    target_basis_d2,
    source_pairing_d1_over_3,
    source_pairing_d2,
):
    """Complete Algorithm 7 for one candidate from Steps 1--4 without intermediate validation."""
    EA, PA, QA = public_key
    E1, _, _, E2, _, _ = ct
    images = recover_torsion_images(candidate, ct)

    timing = {
        "psi2_hd": 0.0,
        "psi2_dual_kernel": 0.0,
        "psi1_d1_over_3_hd": 0.0,
        "psi1_dual_kernel": 0.0,
    }

    # Construct the HD representation of psi2.
    t0 = time.time()
    F2, _ = build_direct_embedding(
        theta4,
        PA,
        QA,
        images["psi2_P"],
        images["psi2_Q"],
        ctx.d2,
        ctx.b,
        kani_parameters=kani2,
    )
    timing["psi2_hd"] += time.time() - t0

    t0 = time.time()
    P2_target, Q2_target = target_basis_d2
    dual2_P, dual2_Q = eval_dual_two_consistent(
        theta4, F2, EA, E2, P2_target, Q2_target
    )
    t = compute_canonical_kernel(
        dual2_P,
        dual2_Q,
        ctx.d2,
        basis=source_basis_d2,
        ePQ=source_pairing_d2,
    )
    timing["psi2_dual_kernel"] += time.time() - t0

    out = []
    path_failures = []
    q1 = ZZ(psi1_degree_d1_over_3)

    for path in psi1_paths:
        path_started = time.time()
        try:
            psi1_d1_over_3_P, psi1_d1_over_3_Q = psi1_d1_over_3_images(ctx, images, path)

            t0 = time.time()
            F1, _ = build_direct_embedding(
                theta4,
                ctx.Pb,
                ctx.Qb,
                psi1_d1_over_3_P,
                psi1_d1_over_3_Q,
                q1,
                ctx.b,
                kani_parameters=kani1_d1_over_3,
            )
            timing["psi1_d1_over_3_hd"] += time.time() - t0

            t0 = time.time()
            P1_target, Q1_target = path["target_basis"]
            dual1_P, dual1_Q = eval_dual_two_consistent(
                theta4, F1, ctx.E0, path["Epre"], P1_target, Q1_target
            )
            s0 = ZZ(
                compute_canonical_kernel(
                    dual1_P,
                    dual1_Q,
                    q1,
                    basis=source_basis_d1_over_3,
                    ePQ=source_pairing_d1_over_3,
                )
            )
            timing["psi1_dual_kernel"] += time.time() - t0

            # Recover s modulo d1/3 and carry all three lifts to Step 10.
            for lift in range(3):
                s = ZZ(s0 + lift * q1) % ZZ(ctx.d1)
                recovered = recover_message_from_preimage(ctx, s, t, candidate["T"])
                recovered.update(
                    {
                        "delta": candidate["delta"],
                        "alpha": candidate["alpha"],
                        "beta": candidate["beta"],
                        "alpha_fault": candidate["alpha_fault"],
                        "images": images,
                        "sign_key": candidate["sign_key"],
                        "psi1_path_index": path["path_index"],
                        "psi1_lift": lift,
                        "s_mod_d1_over_3": s0,
                        "path_elapsed": time.time() - path_started,
                    }
                )
                out.append(recovered)
        except Exception as exc:
            path_failures.append(
                (path["path_index"], type(exc).__name__, str(exc))
            )

    return out, path_failures, timing


def public_attack(ctx, pk_compressed, ct, fault_bit=1):
    """Run Algorithm 7 for all candidates from Steps 1--4 without intermediate validation."""
    total_started = time.time()
    theta4 = load_theta4()

    t0 = time.time()
    public_key, candidates = recover_secret_scalars(
        ctx, pk_compressed, ct, bit_flip_delta_candidates(fault_bit)
    )
    scalar_elapsed = time.time() - t0

    # Precomputation for the HD representations.
    t0 = time.time()
    EA, _, _ = public_key
    E1, _, _, E2, _, _ = ct

    q1, psi1_paths = prepare_psi1_d1_over_3_paths(ctx, E1)
    kani1_d1_over_3 = find_kani_parameters(q1)
    kani2 = find_kani_parameters(ctx.d2)

    source_basis_d1_over_3 = (3 * ctx.Pd1, 3 * ctx.Qd1)
    source_pairing_d1_over_3 = weil_pairing_pari(
        source_basis_d1_over_3[0], source_basis_d1_over_3[1], q1
    )

    source_basis_d2 = torsion_basis(EA, ctx.d2)
    target_basis_d2 = torsion_basis(E2, ctx.d2, even_power=ctx.b)
    source_pairing_d2 = weil_pairing_pari(
        source_basis_d2[0], source_basis_d2[1], ctx.d2
    )
    precompute_elapsed = time.time() - t0

    recovered_results = []
    failures = []
    aggregate_timing = {
        "psi2_hd": 0.0,
        "psi2_dual_kernel": 0.0,
        "psi1_d1_over_3_hd": 0.0,
        "psi1_dual_kernel": 0.0,
    }

    for i, candidate in enumerate(candidates):
        started = time.time()
        try:
            results, path_failures, candidate_timing = complete_public_candidate(
                ctx,
                theta4,
                public_key,
                ct,
                candidate,
                kani1_d1_over_3,
                kani2,
                q1,
                psi1_paths,
                source_basis_d1_over_3,
                source_basis_d2,
                target_basis_d2,
                source_pairing_d1_over_3,
                source_pairing_d2,
            )
            elapsed = time.time() - started
            for result in results:
                result["candidate_index"] = i
                result["candidate_elapsed"] = elapsed
                recovered_results.append(result)
            for pf in path_failures:
                failures.append((i,) + pf)
            for key, value in candidate_timing.items():
                aggregate_timing[key] += value
        except Exception as exc:
            failures.append((i, None, type(exc).__name__, str(exc)))

    return {
        "public_key": public_key,
        "secret_scalar_candidates": candidates,
        "results": recovered_results,
        "failures": failures,
        "psi1_d1_over_3_paths": len(psi1_paths),
        "timing": {
            "scalar_recovery": scalar_elapsed,
            "public_precomputation": precompute_elapsed,
            **aggregate_timing,
            "total": time.time() - total_started,
        },
    }

def _result_key(result):
    """Canonical result modulo <-I_2>."""
    return (
        ZZ(result["s"]),
        ZZ(result["t"]),
        _matrix_key(result["T"]),
        ZZ(result["message"]),
    )


def validate_correct_guesses(alice, ct, attack_result):
    """Validate the correctness of the guesses using the secret key only at the end.

    The secret key is not used for intermediate validation.
    """
    E1, _, S1, E2, R2, S2 = ct
    correct_guesses = []
    validation_rows = []

    groups = {}
    for result in attack_result["results"]:
        groups.setdefault(result["candidate_index"], []).append(result)

    for candidate_index, results in groups.items():
        representative = results[0]
        honest_R1 = representative["alpha"] * representative["images"]["psi1_P"]
        repaired_ct = (E1, honest_R1, S1, E2, R2, S2)

        try:
            s_sk, t_sk, T_sk = alice._trapdoor_inverse(repaired_ct)
            T_sk = canonical_matrix(T_sk)
            sk_message = recover_message_from_preimage(alice, s_sk, t_sk, T_sk)
        except Exception as exc:
            validation_rows.append(
                {
                    "candidate_index": candidate_index,
                    "matches": 0,
                    "error": "%s: %s" % (type(exc).__name__, exc),
                }
            )
            continue

        matching_results = []
        for result in results:
            matches = (
                ZZ(s_sk) == ZZ(result["s"])
                and ZZ(t_sk) == ZZ(result["t"])
                and _matrix_key(T_sk) == _matrix_key(result["T"])
                and ZZ(sk_message["message"]) == ZZ(result["message"])
            )
            if matches:
                matching_results.append(result)
                correct_guesses.append(result)

        validation_rows.append(
            {
                "candidate_index": candidate_index,
                "matches": len(matching_results),
                "s_sk": ZZ(s_sk),
                "t_sk": ZZ(t_sk),
                "T_sk": T_sk,
                "message_sk": ZZ(sk_message["message"]),
            }
        )

    correct_solutions = {_result_key(r) for r in correct_guesses}
    return {
        "rows": validation_rows,
        "matches": correct_guesses,
        "correct_solutions": correct_solutions,
        "unique": len(correct_solutions) == 1,
    }

def run_attack(parameter_set="128", seed=1, fault_bit=1):
    if parameter_set not in PARAMETERS:
        raise ValueError("parameter_set must be one of %s" % sorted(PARAMETERS))

    total_started = time.time()
    set_random_seed(seed)
    random.seed(seed)
    params = PARAMETERS[parameter_set]

    print("[1/5] Initialising FESTA-%s and generating receiver key..." % parameter_set)
    t0 = time.time()
    alice = FESTA(params)
    bob = FESTA(params)
    alice.keygen()
    pk = alice.export_public_key()
    keygen_elapsed = time.time() - t0

    # Test only: the original message is not used by the attack.
    message = ZZ(random.randrange(0, 1 << int(params["lambda_security"])))

    print("[2/5] Producing one ciphertext with a bit-%d fault in alpha..." % fault_bit)
    t0 = time.time()
    ct = bob.encrypt(pk, message, fault_coeff="a00", fault_bit=fault_bit)
    encryption_elapsed = time.time() - t0

    print("[3/5] Running public-only masking/HD/kernel/message recovery...")
    started = time.time()
    attack_result = public_attack(bob, pk, ct, fault_bit=fault_bit)
    public_elapsed = time.time() - started

    print("      secret-scalar candidates: %d" % len(attack_result["secret_scalar_candidates"]))
    print("      reverse psi1 paths per candidate: %d" % attack_result["psi1_d1_over_3_paths"])
    print("      path/lift candidates: %d" % len(attack_result["results"]))
    print("      arithmetic/HD failures: %d" % len(attack_result["failures"]))
    for r in attack_result["results"]:
        print(
            "      candidate %d/path %d/lift %d: delta=%s, padding=%s, m=%s"
            % (
                r["candidate_index"],
                r["psi1_path_index"],
                r["psi1_lift"],
                r["delta"],
                r["padding_ok"],
                r["message"],
            )
        )
    print("      timing breakdown:")
    for key, value in attack_result["timing"].items():
        print("        %-24s %.3f s" % (key + ":", value))

    print("[4/5] Validating correct guesses with the secret key...")
    t0 = time.time()
    validation = validate_correct_guesses(alice, ct, attack_result)
    validation_elapsed = time.time() - t0
    print("      secret-key matches: %d" % len(validation["matches"]))
    print("      distinct correct solutions: %d" % len(validation["correct_solutions"]))
    print("      unique modulo +/- masking: %s" % validation["unique"])

    print("[5/5] Post-Attack Test-validity check")
    recovered_messages = {ZZ(r["message"]) for r in validation["matches"]}
    print("      original message: %s" % message)
    print("      matching recovered messages: %s" % sorted(recovered_messages))
    print("      benchmark summary:")
    print("        receiver keygen:          %.3f s" % keygen_elapsed)
    print("        faulted encryption:       %.3f s" % encryption_elapsed)
    print("        public attack:            %.3f s" % public_elapsed)
    print("        validation of correct guesses: %.3f s" % validation_elapsed)
    print("        full experiment:          %.3f s" % (time.time() - total_started))

    if not validation["unique"]:
        raise AssertionError("validation did not yield a unique correct solution")
    if recovered_messages != {message}:
        raise AssertionError("unique secret-key-validated recovery did not match plaintext")

    return attack_result, validation


if __name__ == "__main__":
    parameter_set = sys.argv[1] if len(sys.argv) > 1 else "128"
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    fault_bit = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    run_attack(parameter_set=parameter_set, seed=seed, fault_bit=fault_bit)