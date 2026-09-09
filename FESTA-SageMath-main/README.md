# FESTA fault-attack proof of concept

This directory contains the FESTA message-recovery fault attack accompanying the paper 
*Fault Injection Attacks on Torsion Masking*. 
The PoC is based on the reference FESTA SageMath implementation.

The attack driver is `faulty_festa.py`. It follows the same dimension-4 Kani/HD path used by the POKE Level-5 attack in this artifact and reuses the bundled `POKE-Level5-Fault-Attack/Theta_dim4_sage` implementation.

## Fault model used by the PoC

The default experiment flips **bit 1** of the diagonal masking coefficient `alpha` in the computation of `R1` only. Since FESTA samples `alpha` as an odd unit modulo `2^b`, flipping bit 1 gives `Delta = +/-2` and keeps both `alpha` and the faulted `alpha'` invertible. The attack solves the resulting linear congruence modulo `2^b`; it does not divide by the even `Delta`.

## Attack phases
The underlying implementation of FESTA is due to https://github.com/FESTA-PKE/FESTA-SageMath
1. 
  **Fault Injection Simulation** From the public key and one faulted ciphertext, the code recovers all masking candidates and unmasks the `2^b`-torsion action of `psi1` and `psi2`. Since FESTA-128 has `d1 = 1 (mod 4)`, a direct `KaniEndoHalf` representation of `psi1` is impossible with the two-square norm equation used by the dimension-4 code. 
  The PoC therefore enumerates the four public reverse 3-isogenies from `E1`, represents the degree-`d1/3` prefix, recovers `s mod d1/3`, and carries all three lifts modulo `d1` to the end. `psi2` is represented directly. 
  It then evaluates the dual isogenies, applies FESTA's `compute_canonical_kernel`, and recovers
  `(s,t,T)` and the plaintext. 
  The receiver secret key, honest encryption trace, and original plaintext are not available to these functions.
2. 
   **Final secret-key Check.** Only after every public candidate has been computed through the complete attack, the script uses the real receiver secret key. 
   Each candidate predicts the honest pre-fault `R1`; the script repairs that one point and invokes FESTA's real trapdoor inverse. 
   This checks which completed recovery is the real one and whether it is unique modulo FESTA's global `+/-` masking equivalence.


## Usage

This PoC currently supports only FESTA-128 end to end. 
From this directory, run:

```text
sage --python -O faulty_festa.py 128 1 1
```
(faulty_festa.py <parameter-set> <seed> <fault-bit>)

The ciphertext is intentionally kept as uncompressed curve/point objects in this PoC so that the injected physical fault is represented directly on the published point. 
The ordinary FESTA `decrypt` method accepts both these uncompressed tuples and the reference compressed representation.

## Example output

A successful FESTA-128 run with seed `1` produces output of the following form.
Exact timings depend on the machine.

```
[1/5] Initialising FESTA-128 and generating receiver key...
[2/5] Producing one ciphertext with a bit-1 fault in alpha...
[3/5] Running public-only masking/HD/kernel/message recovery...
   secret-scalar candidates: 
   reverse psi1 paths per candidate: 
   path/lift candidates: 
   arithmetic/HD failures: 
   candidate 0/path 1/lift 0: delta=-2, padding=False, m=
   candidate 0/path 1/lift 1: delta=-2, padding=False, m=
   candidate 0/path 1/lift 2: delta=-2, padding=False, m=
   candidate 1/path 1/lift 0: delta=-2, padding=False, m=
   candidate 1/path 1/lift 1: delta=-2, padding=False, m=
   candidate 1/path 1/lift 2: delta=-2, padding=True, m=
   timing breakdown:
    scalar_recovery:     s
    public_precomputation:  s
    psi2_hd:         s
    psi2_dual_kernel:    s
    psi1_d1_over_3_hd:    s
    psi1_dual_kernel:    s
    total:          s
**Remark:**
  We remark that, while the recovered torsion information is sufficient for an HD representation of \(\psi_1\), the Kani condition used in our implementation cannot be satisfied directly for \(d_1\equiv1\pmod4\). Hence, we first take one reverse \(3\)-isogeny step from \(E_1\) and work with the resulting degree-\(d_1/3\) isogeny.
[4/5] Validating correct guesses with the secret key...
   secret-key matches: 1
   distinct correct solutions: 1
   unique modulo +/- masking: True
[5/5] Test-validity check (not used by the attack)...
   original message: 151557408999110657826917604970069258585
   matching recovered messages: [151557408999110657826917604970069258585]
   benchmark summary:
    receiver keygen:     s
    faulted encryption:    s
    public attack:      s
    validation of correct guesses: s
    full experiment:     s
```

The reported arithmetic/HD failures are discarded candidate paths, not a
failure of the attack. Success is indicated by a unique matching recovery
class and equality between the original and recovered messages.

## FESTA-128 HD detail

The dimension-4 Kani helper requires `q + a1^2 + a2^2 = 2^e`. For the FESTA-128 `d1`, `d1 = 1 (mod 4)`, so this equation has no solution modulo 4. 
The attack therefore strips one public degree-3 factor on the codomain side.
The resulting prefix degree `q1=d1/3` is `3 (mod 4)` and is compatible with the bundled `KaniEndoHalf` implementation. 
No reverse path or lift is validated during the public phase; all candidates are carried forward unless the requested arithmetic construction itself is undefined.



