#!/usr/bin/env python3
import os, re, sys, math, datetime
import numpy as np
import socket as st

mu_feedback = [0.02]
mu_stress_cue_influx_combis = [[0.02,0.0]]
mu_influx = [0.02]
sdmu = [ 0.02]

# switch rate from P to NP
#s_P2NP = [[0.01,0.02],[ 0.04,0.08], [ 0.1,0.2]]
#s_NP2P = [[0.005,0.01], [ 0.02,0.04], [ 0.05, 0.1]]

s_P2NP_pre = [ 0.1, 0.5, 0.9 ]
s_NP2P_pre = [ 0.1, 0.5, 0.9 ]

s_P2NP = []
s_NP2P = []

for s1_i in s_P2NP_pre:
    s_P2NP += [[s1_i,s1_i]]

for s1_i in s_NP2P_pre:
    s_NP2P += [[s1_i,s1_i]]

#cue_P = [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]
#cue_NP = [ 0.0, 0.2, 0.3, 0.4 ]

cue_P = [ 0.0 ]
cue_NP = [ 0.0 ]

# power of the damage cost function
ad = [ 0.5, 1.0, 1.5 ]

# power of the hormone level survival function
aP = [ 1.0 ]

damage_decay = [ 1.0 ]
damage_due_to_hormone = [ 1.0 ]

# number of replicates
nrep = 3

ctr = 0

init_feedback = [0]
init_stress_influx = [0]
init_stress_baseline_influx = [0]
init_cue_influx = [0]
init_influx = [0]

# attack probability 
p_att = [ 0.5 ]
pleiotropy = [ 0.0 ]
maxtime = "200000"

# attack probability
mort_background = 0.001

exe = "./xstress"


host = st.gethostname()

background = True 
if re.search("anthoxanthum",host) is not None and background:
    background = False

background_str = " & " if background else ""

date = datetime.datetime.now()
base_name = "sim_stress_" + f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"

# make all permutations of parameter combinations
for rep_i in range(0,nrep):
    for mu_feedback_i in mu_feedback:
        for stress_cue_influx_i in mu_stress_cue_influx_combis:
    
            mu_stress_influx_i = stress_cue_influx_i[0]
            mu_cue_influx_i = stress_cue_influx_i[1]

            for mu_influx_i in mu_influx:
                for sdmu_i in sdmu:
                    for s_P2NP_i in s_P2NP:
                        for s_NP2P_i in s_NP2P:
                            for cue_P_i in cue_P:
                                for cue_NP_i in cue_NP:
                                    for p_att_i in p_att:
                                        for ad_i in ad:
                                            for aP_i in aP:
                                                for r_i in damage_decay:
                                                    for u_i in damage_due_to_hormone:
                                                        for init_feedback_i in init_feedback:
                                                            for init_stress_influx_i in init_stress_influx:
                                                                for init_cue_influx_i in init_cue_influx:
                                                                    for init_influx_i in init_influx:
                                                                        for init_stress_baseline_influx_i in init_stress_baseline_influx:
                                                                            for pleiotropy_i in pleiotropy:

                                                                                ctr += 1
                                                                                print("echo " + str(ctr))

                                                                                base_name_i = base_name + "_"  + str(ctr)
                                                                                
                                                                                print(exe + " " 
                                                                                        + str(mu_feedback_i) + " " 
                                                                                        + str(mu_cue_influx_i) + " " 

                                                                                        + str(mu_stress_influx_i) + " " 
                                                                                        + str(mu_influx_i) + " " 

                                                                                        + str(sdmu_i) + " " 
                                                                                        + str(s_P2NP_i[0]) + " " 

                                                                                        + str(s_P2NP_i[1]) + " " 
                                                                                        + str(init_feedback_i) + " " 

                                                                                        + str(init_cue_influx_i) + " " 
                                                                                        + str(init_stress_influx_i) + " " 
                                                                                        + str(init_stress_baseline_influx_i) + " " 
                                                                                        + str(init_influx_i) + " " 
                                                                                        + str(cue_P_i) + " " 

                                                                                        + str(cue_NP_i) + " " 
                                                                                        + str(ad_i) + " " 

                                                                                        + str(aP_i) + " " 
                                                                                        + str(r_i) + " " 

                                                                                        + str(u_i) + " " 
                                                                                        + str(mort_background) + " " 

                                                                                        + str(p_att_i) + " " 
                                                                                        + str(pleiotropy_i) + " " 
                                                                                        + str(maxtime) + " "
                                                                                        + str(0) + " " 

                                                                                        + base_name_i + " "
                                                                                        + str(background_str)
                                                                                        )
