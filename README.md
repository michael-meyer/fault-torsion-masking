# Fault Injection Attacks on Torsion Masking

This repository contains the proof-of-concept simulations for the fault injection attacks presented in *Fault Injection Attacks on Torsion Masking*.

The artifact contains implementations of the attacks on:

- PRISM;
- POKÉ; and
- FESTA.

The simulations are based on the reference SageMath implementations of the corresponding schemes. We modify the torsion masking computations to simulate the fault models described in the paper and implement the corresponding post-processing for key or message recovery.

Each directory contains the code and instructions for the corresponding attack. Please refer to the scheme-specific `README.md` files for implementation details and commands for reproducing the experiments.