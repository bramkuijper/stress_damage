#!/usr/bin/env python3
import numpy as np
import os.path

autocorr = [ 0, 0.1, 0.3, 0.5, 0.7, 0.9 ]
risk = [ 0.05, 0.1, 0.2 ]

# probabilities of leaving and arriving

pLA = []

for autocorr_i in autocorr:
    for risk_i in risk:
        pL = (1.0 - autocorr_i)/(1.0+(risk_i/(1.0-risk_i)));
        pA = 1.0 - pL - autocorr_i

        pLA += [[round(pL,4),round(pA,4)]]

Kfec = [ 0.05 ]
Kmort = [ 0.01 ]

pAttack = [0.5]

alpha = 1.0


exe = "stress_damage_lh.exe"

suppress_echo = True 

the_dir = "./"

background = False

ctr = 1

for pLA_i in pLA:
    for Kfec_i in Kfec:
        for Kmort_i in Kmort:
            for pAttack_i in pAttack:

                if not suppress_echo:
                    print("echo " + str(ctr))

                print(os.path.join(the_dir,exe) + " " +
                        str(pLA_i[0]) + " " +
                        str(pLA_i[1]) + " " +
                        str(pAttack_i) + " " +
                        str(alpha) + " " +
                        str(Kmort_i) + " " +
                        str(Kfec_i) + " " +
                        ("&" if background else "")
                        )
                ctr+=1
