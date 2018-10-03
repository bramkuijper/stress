#!/usr/bin/env python3
import os, re, sys, math
from numpy import *

mu_feedback = [0, 0.01];
mu_stress_influx = [0, 0.01];
mu_influx = [0, 0.01];
sdmu = [ 0.01 ];

# switch rate from P to NP
s_P2NP = [[ 0.1 ,0.5]]
s_NP2P = [[ 0.01, 0.25 ]]

s_12 = [[ 0.1,0.1] ]

cue_P = [ 0.8 ]
cue_NP = [ 0.08 ]

s0 = [ 0.5 ]
ad = [ 0.5, 1.0, 2.0 ]
aP = [ 0.5, 1.0, 2.0 ]

dmax = [ 100 ]
zmax = [ 100 ]

r = [ 0.1, 0.5, 0.9 ]
u = [ 0.5, 1, 10 ]

# number of replicates
nrep = 3

ctr = 0

init_feedback = 0.0
init_stress_influx = 0.0
init_influx = 0.0

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
                                                            for r_i in r:
                                                                for u_i in u:

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
