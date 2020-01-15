#!/usr/bin/env python3
import os, re, sys, math
import numpy as np
import socket as st

mu_feedback = [0.02]
mu_stress_influx = [0.02]
mu_influx = [0.02]
sdmu = [ 0.02]

# switch rate from P to NP
#s_P2NP = [[0.01,0.02],[ 0.04,0.08], [ 0.1,0.2]]
#s_NP2P = [[0.005,0.01], [ 0.02,0.04], [ 0.05, 0.1]]

s_P2NP_pre = [ 0.01 ]
s_NP2P_pre = [ 0.01 ]

s_P2NP = []
s_NP2P = []

for s1_i in s_P2NP_pre:
    s_P2NP += [[s1_i,s1_i]]

for s1_i in s_NP2P_pre:
    s_NP2P += [[s1_i,s1_i]]

# from one world to another
s_12 = [[ 0.1,0.1]]

#cue_P = [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]
#cue_NP = [ 0.0, 0.2, 0.3, 0.4 ]

cue_P = [ 0.0 ]
cue_NP = [ 0.0 ]

s0 = [ 0.1]
ad = [ 2 ]
aP = [ 0.5 ]

dmax = [ 1 ]
zmax = [ 1 ]

damage_decay = [ 1.0 ]

damage_due_to_hormone = [ 0.75 ]

# number of replicates
nrep = 20

ctr = 0

init_feedback = [0]
init_stress_influx = [0]
init_influx = [0]

mort_background = 0.1

exe = "./xstress"


host = st.gethostname()

background = False
if re.search("anthoxanthum",host) is not None:
    background = True

background_str = " & " if background else ""

# make all permutations of parameter combinations
for rep_i in range(0,nrep):
    for mu_feedback_i in mu_feedback:
        for mu_stress_influx_i in mu_stress_influx:
            for mu_influx_i in mu_influx:
                for sdmu_i in sdmu:
                    for s_P2NP_i in s_P2NP:
                        for s_NP2P_i in s_NP2P:
                            for s_12_i in s_12:
                                for cue_P_i in cue_P:
                                    for cue_NP_i in cue_NP:
                                        for s0_i in s0:
                                            for ad_i in ad:
                                                for aP_i in aP:
                                                    for dmax_i in dmax:
                                                        for zmax_i in zmax:
                                                            for r_i in damage_decay:
                                                                for u_i in damage_due_to_hormone:
                                                                    for init_feedback_i in init_feedback:
                                                                        for init_stress_influx_i in init_stress_influx:
                                                                            for init_influx_i in init_influx:

                                                                                ctr += 1
                                                                                print("echo " + str(ctr))

                                                                                print(exe + " " 
                                                                                        + str(mu_feedback_i) + " " 
                                                                                        + str(mu_stress_influx_i) + " " 
                                                                                        + str(mu_influx_i) + " " 
                                                                                        + str(sdmu_i) + " " 
                                                                                        + str(s_P2NP_i[0]) + " " 
                                                                                        + str(s_P2NP_i[1]) + " " 
                                                                                        + str(s_NP2P_i[0]) + " " 
                                                                                        + str(s_NP2P_i[1]) + " " 
                                                                                        + str(s_12_i[0]) + " " 
                                                                                        + str(s_12_i[1]) + " " 
                                                                                        + str(init_feedback_i) + " " 
                                                                                        + str(init_stress_influx_i) + " " 
                                                                                        + str(init_influx_i) + " " 
                                                                                        + str(cue_P_i) + " " 
                                                                                        + str(cue_NP_i) + " " 
                                                                                        + str(s0_i) + " " 
                                                                                        + str(ad_i) + " " 
                                                                                        + str(aP_i) + " " 
                                                                                        + str(dmax_i) + " " 
                                                                                        + str(zmax_i) + " " 
                                                                                        + str(r_i) + " " 
                                                                                        + str(u_i) + " " 
                                                                                        + str(mort_background) + " " 
                                                                                        + str(background_str)
                                                                                        )
