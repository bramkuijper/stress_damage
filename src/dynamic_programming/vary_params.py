#!/usr/bin/env python3
import numpy as np
import os.path, sys
import datetime


# this script is used to vary command-line parameters 
# that will be fed into the simulation

# specify values of the autocorrelation and risk

#autocorr = np.linspace(start=0,stop=0.9,num=50)
#risk = [ 0.02,0.05,0.1,0.2,0.3,0.5 ]

autocorr = [ 0.7,0.9 ]
risk = [ 0.05 ]

n_repro_bout = [1]

# translate those in probabilities
# of leaving and arriving
pLA = []

for autocorr_i in autocorr:
    for risk_i in risk:
        pL = (1.0 - autocorr_i)/(1.0+(risk_i/(1.0-risk_i)));
        pA = 1.0 - pL - autocorr_i

        pLA += [[round(pL,4),round(pA,4)]]

Kfec = [ 0.1 ]
Kmort = [ 0.00 ]

damage_terminal_rewards = [0,1]

pAttack = [0.5]

alpha = 1.0


repair = [0.0]
#repair = np.linspace(start=0,stop=10,num=50)

exe = "stress_damage_terminal.exe"
#exe = "stress_damage_lh.exe"
#exe = "stress_damage.exe"

suppress_echo = True 

the_dir = "./"

background = False

ctr = 1

full_exe_name = os.path.join(the_dir,exe)

nrep = 1

# standard exe
if exe == "stress_damage.exe":
    for autocorr_i in autocorr:
        for risk_i in risk:
            for repair_i in repair:
                for rep_i in range(0,nrep):
                    for pattack_i in pAttack:
                        print(
                                f"{full_exe_name} {autocorr_i} "
                                + f"{risk_i} {repair_i} "
                                + f"{rep_i} "
                                + f"{pattack_i} "
                                )

                        ctr+=1

    # exiting the non-lh loop
    sys.exit(1)

date = datetime.datetime.now()
base_name = "output_file" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"


if exe == "stress_damage_terminal.exe":
    for autocorr_i in autocorr:
        for risk_i in risk:
            for repair_i in repair:
                for rep_i in range(0,nrep):
                    for pattack_i in pAttack:
                        for Kfec_i in Kfec:
                            for damage_terminal_rewards_i in damage_terminal_rewards:
                                print(
                                        f"{full_exe_name} " 
                                        + f"{autocorr_i} "
                                        + f"{risk_i} "
                                        + f"{repair_i} "
                                        + f"{rep_i} "
                                        + f"{pattack_i} "
                                        + f"{Kfec_i} "
                                        + f"{damage_terminal_rewards_i} "
                                        + base_name + "_" + str(ctr) 
                                        )
                                ctr+=1

    # exiting the non-lh loop
    sys.exit(1)
    
# finally the same for the life-history data:
for pLA_i in pLA:
    for Kfec_i in Kfec:
        for Kmort_i in Kmort:
            for pAttack_i in pAttack:
                for repair_i in repair:
                    for repro_bout_i in n_repro_bout:
                        if not suppress_echo:
                            print("echo " + str(ctr))

                        print(full_exe_name + " " +
                                str(pLA_i[0]) + " " +
                                str(pLA_i[1]) + " " +
                                str(pAttack_i) + " " +
                                str(alpha) + " " +
                                str(Kmort_i) + " " +
                                str(Kfec_i) + " " +
                                str(repair_i) + " " +
                                str(repro_bout_i) + " " +
                                base_name + "_" + str(ctr) +
                                ("&" if background else "")
                                )
                        ctr+=1
