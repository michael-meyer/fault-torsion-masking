from prism import PRISM, PRISM_verify
from sage.all import (
    pari, Zmod, matrix, ZZ, isqrt, gcd, Integers, floor, discrete_log
)
import random


try:
    from . import params, ec, hd, misc
    from . import quaternions as qt
    from . import qlapoti as qlpt
except ImportError:
    import params, ec, hd, misc
    import quaternions as qt
    import qlapoti as qlpt

lvl = 1
params.set_prism_params(lvl)


'''
specify the fault type below

allowed types:     'bit_flip' flips last bit
                   'random_byte' assigns a random value to the last byte (the simulation
                        only faults the last three bits instead of the full byte to keep
                        runtimes reasonable in this proof of concept)
targetting a coefficient: allowed options '11', '12', '21', '22' to target m_ij
'''

#fault = 'bit_flip'
fault = 'random_byte'

#candidates for delta with targeted bit-flip of last bit: 1 and -1
if fault == 'bit_flip':
    delta1_candidates = [-1, 1]
    delta2_candidates = [-1, 1]

if fault == 'random_byte':
    
    # RANDOM LAST BYTE
    #delta1_candidates = [2*k+1 for k in range(-128,0)]
    #delta2_candidates = [2*k+1 for k in range(-128,128)]
    
    # RANDOM LAST 3 BITS
    delta1_candidates = [2*k+1 for k in range(-8,0)]
    delta2_candidates = [2*k+1 for k in range(-8,8)]


# run keygen
alice = PRISM()

#fault-less run (only used to double-check results)
msg = 'Hello world'
sigma = alice.sign(msg)

#make sure verifiaction works for a valid signature
assert PRISM_verify(msg, sigma, alice.pk)

#secret masking scalars and isogeny degree; used only to double-check results of our attack at the end
sk_m11, sk_m12, sk_m21, sk_m22 = alice.sk[-1].list()
sk_deg = (alice.sk[-2]).norm()

#public key points
Pvk = alice.pk[1][0]
Qvk = alice.pk[1][1]



if __name__ == "__main__":

    #counter for counting faulty runs
    #since the simulation repeats faulty runs if delta is even
    num_faults = 0

    #filter out even delta values
    delta1_odd = False
    delta2_odd = False

    while delta1_odd == False or delta2_odd == False: 

        #first faulty run to get P'
        if delta1_odd == False:
            sigma_fault1 = alice.sign(msg, fault, '11')
            num_faults += 1

        #second faulty run to get Q'
        if delta2_odd == False:
            sigma_fault2 = alice.sign(msg, fault, '22')
            num_faults += 1


        #get correct points and assert that faults were working
        P_fault = sigma_fault1[1][0]
        Q_fault = sigma_fault2[1][1]
        P = sigma_fault2[1][0]
        Q = sigma_fault1[1][1]

        P_diff = P - P_fault
        Q_diff = Q - Q_fault

        #check point order to determine if delta is odd
        #if not, we repeat faulty runs
        if P_diff.order() == 2**params.a:
            delta1_odd = True
        else:
            print("P_diff does not have full order, trying a new faulty run")
        if Q_diff.order() == 2**params.a:
            delta2_odd = True
        else:
            print("Q_diff does not have full order, trying a new faulty run")
            
    assert P == sigma[1][0]
    assert Q == sigma[1][1]
    assert P != P_fault
    assert Q != Q_fault
    assert(P.order() == 2**params.a)
    assert(Q.order() == 2**params.a)

    #loop through possible combinations of delta1 and delta2
    #following Algorithm 2

    found = False 

    for d1 in delta1_candidates:

        if found == True:
            break

        print(f"checking for all combinations with delta1 = {d1}")

        for d2 in delta2_candidates:

            if found == True:
                break


            delta1 = Zmod(2**params.a)(d1)
            delta2 = Zmod(2**params.a)(d2)

            #compute inverses for multiplication in Step 1
            delta1_inv = delta1 ** (-1)
            delta2_inv = delta2 ** (-1)


            #compute phi(P0) and phi(Q0) from Step 1
            phi_P0 = delta1_inv * P_diff
            phi_Q0 = delta2_inv * Q_diff
            assert(phi_Q0.order() == 2**params.a)
            assert(phi_P0.order() == 2**params.a)
            

            #compute pairings from Step 2
            lam_PphiP0 = P.weil_pairing(phi_P0, 2**params.a)
            lam_PphiQ0 = P.weil_pairing(phi_Q0, 2**params.a)
            lam_QphiP0 = Q.weil_pairing(phi_P0, 2**params.a)
            lam_QphiQ0 = Q.weil_pairing(phi_Q0, 2**params.a)
            lam_phiP0phiQ0 = phi_P0.weil_pairing(phi_Q0, 2**params.a)
            lam_phiQ0phiP0 = lam_phiP0phiQ0**(-1)


            #compute discrete logs from Step 3
            m11 = discrete_log(lam_PphiQ0, lam_phiP0phiQ0, 2**params.a)
            assert(lam_PphiQ0 == lam_phiP0phiQ0**m11)
            m12 = discrete_log(lam_PphiP0, lam_phiQ0phiP0, 2**params.a)
            assert(lam_PphiP0 == lam_phiQ0phiP0**m12)
            m21 = discrete_log(lam_QphiQ0, lam_phiP0phiQ0, 2**params.a)
            assert(lam_QphiQ0 == lam_phiP0phiQ0**m21)
            m22 = discrete_log(lam_QphiP0, lam_phiQ0phiP0, 2**params.a)
            assert(lam_QphiP0 == lam_phiQ0phiP0**m22)
            

            #recovery of HD representation

            #compute points from Step 4
            M = matrix(Zmod(2**params.a), 2, [m11, m12, m21, m22])
            M_inv = M.inverse()
            m11_inv, m12_inv, m21_inv, m22_inv = M_inv.list()
            phisk_P0 = m11_inv * Pvk + m12_inv * Qvk
            phisk_Q0 = m21_inv * Pvk + m22_inv * Qvk
            assert(phisk_P0.order() == 2**params.a)
            assert(phisk_Q0.order() == 2**params.a)

            #compute pairings from Step 5
            lam_phiskP0phiskQ0 = phisk_P0.weil_pairing(phisk_Q0, 2**params.a)
            lam_P0Q0 = (params.P0).weil_pairing(params.Q0, 2**params.a)

            deg = discrete_log(lam_phiskP0phiskQ0, lam_P0Q0, 2**params.a)


            #check HD representation (see Sec. 4.2, "Validating correct guesses")

            #generate isogeny psi: E0 -> E1 of degree deg2 = 2^a-deg

            deg2 = Zmod(2**params.a)(2**params.a-deg)

            #points to push through psi
            P00 = params.P0
            Q00 = params.Q0

            #generate ideal of norm deg2
            I_psi =  qt.RandomIdealGivenNorm(int(deg2), False)

            #compute corresponding isogeny and evaluate P00, Q00
            E1, P1, Q1 = qlpt.IdealToIsogeny(I_psi,P00,Q00)
            assert(P1.order() == 2**params.a)
            assert(Q1.order() == 2**params.a)

            #compute kernel points
            K = (P1, phisk_P0), (Q1, phisk_Q0)

            try:
                isog = hd.Dim2Iso(K, params.a)
            except ValueError:
                continue

            #check if E0 is in codomain to validate the guess
            if (isog.codomain()[0]).j_invariant() == 1728 or (isog.codomain()[1]).j_invariant() == 1728:
                print("found a product containing E0")
                print(f"delta values: delta1 = {d1}, delta2 = {d2}")
                print("check against secret key values:")
                print("checking m_ij coefficients vs. secret key coefficients:")
                correctness_check = (m11 == sk_m11 or m11 == -sk_m11) and (m12 == sk_m12 or m12 == -sk_m12) and (m21 == sk_m21 or m21 == -sk_m21) and (m22 == sk_m22 or m22 == -sk_m22)
                print(correctness_check)
                # occaisonally there exists two candidates corresponding to solutions with E0 in the codomain
                # in practice all candidate delta1, delta2 solutions should be checked to determine when this occurs
                if correctness_check:
                    print("checking secret isogeny degree:")
                    print(deg == sk_deg)
                    found = True
                    print("required faults:")
                    print(num_faults)


        

        


