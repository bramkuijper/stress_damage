#!/usr/bin/env python3
import numpy as np
import os.path, sys

# specify values of the autocorrelation and risk
#autocorr = [ 0, 0.1, 0.3, 0.5, 0.7, 0.9 ]
#risk = [ 0.05, 0.1, 0.2 ]

autocorr = [ 0.3 ]
risk = [ 0.2 ]

# translate those in probabilities
# of leaving and arriving
pLA = []

for autocorr_i in autocorr:
    for risk_i in risk:
        pL = (1.0 - autocorr_i)/(1.0+(risk_i/(1.0-risk_i)));
        pA = 1.0 - pL - autocorr_i

        pLA += [[round(pL,4),round(pA,4)]]

Kfec = [ 0.01 ]
Kmort = [ 0.001 ]

pAttack = [0.5]

alpha = 1.0


repair = [0.0,1.0]
#repair = np.linspace(start=0,stop=10,num=50)

#exe = "stress_damage_lh.exe"
exe = "stress_damage.exe"

suppress_echo = True 

the_dir = "./"

background = False

ctr = 1

full_exe_name = os.path.join(the_dir,exe)

# standard exe
if "lh" not in exe:
    for autocorr_i in autocorr:
        for risk_i in risk:
            for repair_i in repair:
                print(f"{full_exe_name} {autocorr_i} "
                        + f"{risk_i} {repair_i}")

                ctr+=1
    sys.exit(1)




for pLA_i in pLA:
    for Kfec_i in Kfec:
        for Kmort_i in Kmort:
            for pAttack_i in pAttack:

                if not suppress_echo:
                    print("echo " + str(ctr))

                print(full_exe_name + " " +
                        str(pLA_i[0]) + " " +
                        str(pLA_i[1]) + " " +
                        str(pAttack_i) + " " +
                        str(alpha) + " " +
                        str(Kmort_i) + " " +
                        str(Kfec_i) + " " +
                        ("&" if background else "")
                        )
                ctr+=1
