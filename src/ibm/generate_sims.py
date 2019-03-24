#!/usr/bin/env python3
import os, re, sys, math
import numpy as np

mu_feedback = [0.02]
mu_stress_influx = [0.02]
mu_influx = [0.02]
sdmu = [ 0.02]

# switch rate from P to NP
#s_P2NP = [[0.01,0.02],[ 0.04,0.08], [ 0.1,0.2]]
#s_NP2P = [[0.005,0.01], [ 0.02,0.04], [ 0.05, 0.1]]

s_P2NP_pre = [ 0.01, 0.02, 0.05, 0.075, 0.1]
s_NP2P_pre = [ 0.01, 0.02, 0.05, 0.075, 0.1]

#s_P2NP_pre = [ 0.1 ]
#s_NP2P_pre = [ 0.05 ]


s_P2NP = []
s_NP2P = []

for s1_i in s_P2NP_pre:
    s_P2NP += [[s1_i,s1_i]]

for s1_i in s_NP2P_pre:
    s_NP2P += [[s1_i,s1_i]]


s_12 = [[ 0.1,0.1]]

cue_P = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ]
cue_NP = [ 0.0, 0.05, 0.1, 0.2, 0.3, 0.4 ]

#cue_P = [ 0.8 ]
#cue_NP = [ 0.4 ]

s0 = [ 0.5 ]
ad = [ 0.5 ]
aP = [ 0.5 ]

dmax = [ 100 ]
zmax = [ 100 ]

damage_decay = [ 0.9 ]
damage_due_to_hormone = [ 1.0 ]

# number of replicates
nrep = 1

ctr = 0

init_feedback = 1.0
init_stress_influx = 10.0
init_influx = 2.0

exe = "./xstress"

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
                                                                            + str(init_feedback) + " " 
                                                                            + str(init_stress_influx) + " " 
                                                                            + str(init_influx) + " " 
                                                                            + str(cue_P_i) + " " 
                                                                            + str(cue_NP_i) + " " 
                                                                            + str(s0_i) + " " 
                                                                            + str(ad_i) + " " 
                                                                            + str(aP_i) + " " 
                                                                            + str(dmax_i) + " " 
                                                                            + str(zmax_i) + " " 
                                                                            + str(r_i) + " " 
                                                                            + str(u_i) + " " 
                                                                            )
