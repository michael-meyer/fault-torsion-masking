# POKE Level-5 fault attack

This directory contains the POKE message-recovery fault attack accompanying the paper 
*Fault Injection Attacks on Torsion Masking*. 
This package runs one focused experiment against the original POKE Level-5 parameters (`a=254`, `b=324`).  It combines the paper's torsion-masking fault attack with a five-step reverse search in the 3-isogeny graph.

## Requirements

- WSL/Linux
- SageMath 10.0 or newer (tested with SageMath 10.7)
- `Theta_dim4_sage/` in this directory

Activate the tested Conda environment when applicable:

```bash
conda activate sage
```

## Run

```bash
sage lvl5_paper_attack.sage 1
```

The optional argument is the experiment seed.  The final single-core test on the tested machine took `3169.534 s` (about 52.8 minutes).

The program prints each setup phase immediately.  During the long HD search it updates an in-place progress bar with the checked/maximum candidate count, elapsed time, and a worst-case ETA.  A valid candidate is normally found before the bar reaches 100%.

## Attack flow

1. Generate an original Level-5 key and one faulted ciphertext.
2. Recover four masking-scalar candidates from public pairings.
3. Enumerate the `4*3^4 = 324` non-backtracking reverse paths of length five
   from the public ciphertext curve `EB`.
4. For each scalar/path pair, construct the degree-`3^319` HD candidate.
5. Reject invalid candidates when public HD evaluation fails (normally at the
   theta gluing consistency assertion).
6. For a surviving HD candidate, solve a DLP for `r3 mod 3^319`.
7. Enumerate the remaining `3^5 = 243` lifts and retain the unique tail whose
   codomain is the public `EB`.
8. Rebuild both `psi` and `psi'`, require their codomains to equal public `EB`
   and `EAB`, recover the public mixing matrix, and decrypt.

## No secret correctness oracle

Candidate selection does not read the secret key, the encryption scalar `r3`,
an encryption trace, or the original plaintext.  It uses only:

- public parameters and the faulted ciphertext;
- the stated bit-flip fault model;
- public HD evaluation success/failure;
- public `EB` and `EAB` codomain checks.

The experiment encrypts 32 random bytes and discards the original value.  
It prints only the plaintext recovered by the attack.

## Expected output

```text
[1/6] Loading original Level-5 parameters and Theta_dim4...
...
[5/6] Testing at most 1296 public HD candidates...
    [##########################----] ... 1140/1296 ...
      unique public candidate found after 1140 HD attempts
[6/6] Validating public EB/EAB and decrypting...
masking candidates: 4
reverse paths: 324
HD candidates attempted: 1140
public EB/EAB validation: passed
recovered plaintext (hex):
5e55a2487ea1a1ea296cdd6929b6dbef701a371ededc740923583147f3a46ee8
total time: 3169.534 s
```

## Files required at runtime

```text
lvl5_paper_attack.sage      focused attack entry point
faulty_poke_direct.sage     fault, pairing, and HD helper functions
POKE_PKE.sage               original POKE implementation
montgomery_isogenies/       POKE Kummer/isogeny arithmetic
theta_isogenies/            original POKE theta isogenies
theta_structures/           original POKE theta structures
utilities/                  pairing, DLP, torsion, and arithmetic helpers
Theta_dim4_sage/            external dimension-4 Kani implementation
```

`*.sage.py`, `*.pyc`, and `__pycache__/` are generated caches and are not
required.

## Release archive

The delivery ZIP contains this README, the Apache-2.0 license, 
the three Sage source files, the original POKE arithmetic modules, and 
the minimal `Theta_dim4_sage/pkg` module set exercised by the complete Level-5 attack.


The bundled `Theta_dim4_sage/pkg` code is derived from the experimental
Theta_dim4 SageMath library, copyright (c) 2024 Pierrick Dartois.
