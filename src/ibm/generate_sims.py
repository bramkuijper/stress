#!/usr/bin/env python3
import os, re, sys, math
from numpy import *

mu_feedback = [0, 0.01];
mu_stress_influx = [0, 0.01];
mu_influx = [0, 0.01];
sdmu = [ 0.01 ];

# switch rate from P to NP
s_P2NP = [ 0.1, 0.5 ]
s_NP2P = [ 0.01, 0.5 ]

s_12 = [ 0.1 ]

cue_P = [ 0.8 ]
cue_NP = [ 0.08 ]

s0 = [ 0.5 ]
ad = [ 0.5, 1.0, 2.0 ]
aP = [ 0.5, 1.0, 2.0 ]

dmax = [ 

for rep_i in range(0,nrep):
    for d_i in d:
        for combi_i in fecmortcombis:

            ca1 = 1.0
            cm1 = combi_i[1]
            ca2 = 1.0
            cm2 = combi_i[1]

            Fa1 = combi_i[0]
            Fm1 = 1.0
            Fa2 = combi_i[0]
            Fm2 = 1.0

            for p1_init_i in p1_init:
                for p2_init_i in p2_init:
                    for f2 in freq_patch_2:
                        for sbar_i in sbar:
                            for k_i in k:
                                s2 = sqrt(((1.0-f2)/f2) * 10**(2*sbar_i))
                                s1 = 10**(2*sbar_i) / s2

                                ctr += 1
                                print("echo " + str(ctr))

                                print(exe + " " 
                                        + str(s1) + " " 
                                        + str(s2) + " " 
                                        + str(cm1) + " " 
                                        + str(cm2) + " " 
                                        + str(ca1) + " " 
                                        + str(ca2) + " " 
                                        + str(Fm1) + " "  
                                        + str(Fm2) + " "  
                                        + str(Fa1) + " " 
                                        + str(Fa2) + " "  
                                        + str(k_i) + " " 
                                        + " 0.02 0.01 " 
                                        + str(d_i) + " " 
                                        + str(p1_init_i) + " " 
                                        + str(error_prob) + " " 
                                        + str(p1_init_i) + " " 
                                        + str(p2_init_i) + " " 
                                        )
