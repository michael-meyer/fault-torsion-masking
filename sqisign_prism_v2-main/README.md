# SQIsign / PRISM v2 in Sage

Proof-of-concept simulation for fault attacks on PRISM accompanying the paper 'Fault Injection Attacks on Torsion Masking'.

Based on the PRISM implementation in SageMath from
https://github.com/KULeuven-COSIC/sqisign_prism_v2 


The implementation simulates the fault model by faulting masking coefficients during the torsion point masking.
Supported fault types are flipping the last bit, or assigning a random value to the last byte (resp. to its last 3 bits for faster post-precessing times).

The underlying implementation of PRISM is due to
https://github.com/KULeuven-COSIC/sqisign_prism_v2.
We only modified the PRISM implementation to add the simulation of faults and to obtain a deterministic scheme as described in the initial PRISM paper. 

## Usage - SQIsign

The file `faulty_prism.py` runs the simulation and prints 
the progress and the identified secret data.
It can be run with `sage --python -O faulty_prism.py'

By default, it runs with NIST security level 1 and random_byte faults.
These values can be changed in the first lines of `faulty_prism.py`
as detailed there.
Note, the script occaisionally finds 2 candidate solutions, of which 1 is always correct.

