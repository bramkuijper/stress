#!/usr/bin/env python3
import os
import re
import sys
import math
import datetime
import numpy as np
import socket as st
from scipy import optimize


def logistic(x):
    return(1.0 / (1.0 + np.exp(-x)))

def logit(x):
    return(math.log(x/(1.0 - x)))



exe = "./xstress"


###### mutation rates ######

# mutation rate for the clearance trait
mu_clearance = [0.0005]

# combination of mutation rates for 
# 1. the stress influx and 
# 2. the stress influx slope and 
# 3. the cue influx
mu_stress_cue_influx_combis = [[0.0005,0.0, 0.0]]

# mutation rate of the baseline influx
mu_influx = [0.0005]

# mutation rate of the hstart trait
mu_hstart = [0.0005]

# mutational standard deviation
sdmu = 0.10




###### environment ######


s_P2NP = [ 0.9 ]
s_NP2P = [ 0.1 ]

# switch rates as a combination of autocorrelation and risk
# autocorr (rho = 1 - (sP2NP + sNP2P)
#autocorr = np.arange(0,1,0.1)
autocorr = [0,0.3,0.9]

# risk (r) = sNP2P/(sP2NP + sNP2P)
#risk = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 ]
risk = [0.05, 0.1]


# hence:
# s_P2NP = (1-r) ( 1 - rho)
# s_NP2P = r ( 1 - rho)


s_combinations = []

for rho_i in autocorr:
    for r_i in risk:
        s_P2NP = (1 - r_i) * (1 - rho_i)
        s_NP2P = r_i * (1 - rho_i)
        s_combinations.append([round(s_P2NP,3), round(s_NP2P,3)])

#[print([s[0],s[1],s[1]/(s[0] + s[1]),1 - s[0] - s[1]]) for s in s_combinations]
#sys.exit(1)


# probability of receiving a cue in the predator environment
cue_P = [ 0.0 ]
# probability of receiving a cue in the predator-free environment
cue_NP = [ 0.0 ]

# power of the damage cost function
ad = [ 1.5 ]

# power of the hormone level survival function
aP = [ 1.0 ]

# attack probability 
p_att = [ 0.5 ]

s0 = 0.0

dmax = 1.0
zmax = 1.0
min_clearance = [0.01]


# rate at which damage decays over time
# 1.0: full decay, i.e
# no damage lingering after current timestep
# 0.0: no decay whatsoever
damage_decay = [ 1.0 ]
damage_due_to_hormone = [ 1.0 ]

# number of replicates
nrep = 5

# background mortality
mort_background = 0.002

# number of timesteps
maxtime = 4e06
#maxtime = 3
# number of timesteps stress influx lasts
tmax_stress_influx = [ 1 ]

#### initial values for all traits ####

init_clearance = [0.1]
init_stress_influx = [0.4]
init_stress_influx_slope = [0.0]
init_influx = [0.05]
init_hstart = [0.5]



##### determine whether we should run jobs in background #####
background_jobs = False 

host = st.gethostname()

if re.search("anthoxanthum",host) is None and background_jobs:
    background_jobs = False

background_str = " & " if background_jobs else ""








date = datetime.datetime.now()
base_name = "sim_stress_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"
ctr = 0

# make all permutations of parameter combinations
for rep_i in range(0,nrep):
    for mu_clearance_i in mu_clearance:
        for stress_cue_influx_i in mu_stress_cue_influx_combis:
    
            mu_stress_influx_i = stress_cue_influx_i[0]
            mu_stress_influx_slope_i = stress_cue_influx_i[1]
            mu_cue_influx_i = stress_cue_influx_i[2]

            for mu_influx_i in mu_influx:
                
                for mu_hstart_i in mu_hstart:

                    for s_combi_i in s_combinations:
                        s_P2NP_i = s_combi_i[0]
                        s_NP2P_i = s_combi_i[1]

                        for cue_P_i in cue_P:
                            for cue_NP_i in cue_NP:
                                for p_att_i in p_att:

                                    for ad_i in ad:
                                        for aP_i in aP:

                                            for min_clearance_i in min_clearance:
                                                for r_i in damage_decay:
                                                    for u_i in damage_due_to_hormone:

                                                        for init_clearance_i in init_clearance:
                                                            for init_stress_influx_slope_i in init_stress_influx_slope:
                                                                for init_stress_influx_i in init_stress_influx:
                                                                    for init_influx_i in init_influx:
                                                                        for init_hstart_i in init_hstart:
                                                                            for tmax_stress_influx_i in tmax_stress_influx:

                                                                                ctr += 1
                                                                                print("echo " + str(ctr))

                                                                                base_name_i = base_name + "_"  + str(ctr)

                                                                                print(exe + " " 
                                                                                      + str(mu_clearance_i) + " " 
                                                                                      + str(mu_stress_influx_i) + " " 
                                                                                      + str(mu_stress_influx_slope_i) + " " 
                                                                                      + str(mu_influx_i) + " " 
                                                                                      + str(mu_hstart_i) + " " 
                                                                                
                                                                                      + str(sdmu) + " " 
                                                                                      + str(sdmu) + " " 
                                                                                      + str(sdmu) + " " 
                                                                                      + str(sdmu) + " " 
                                                                                      + str(sdmu) + " " 

                                                                                      + str(s_P2NP_i) + " " 
                                                                                      + str(s_NP2P_i) + " " 
                                                                                
                                                                                      + str(init_clearance_i) + " " 
                                                                                      + str(init_stress_influx_i) + " " 
                                                                                      + str(init_stress_influx_slope_i) + " " 
                                                                                      + str(init_influx_i) + " " 
                                                                                      + str(init_hstart_i) + " " 
                                                                                
                                                                                      + str(cue_P_i) + " " 
                                                                                      + str(cue_NP_i) + " " 
                                                                                
                                                                                      + str(s0) + " " 
                                                                                      + str(ad_i) + " " 
                                                                                      + str(aP_i) + " " 
                                                                                
                                                                                      + str(dmax) + " " 
                                                                                      + str(zmax) + " " 
                                                                                      + str(min_clearance_i) + " " 
                                                                                
                                                                                      + str(r_i) + " " 
                                                                                      + str(u_i) + " " 
                                                                                
                                                                                      + str(mort_background) + " " 
                                                                                      + str(p_att_i) + " " 

                                                                                      + str(maxtime) + " "
                                                                                      + str(tmax_stress_influx_i) + " "
                                                                                      + str(0) + " " 
                                                                                
                                                                                      + base_name_i + " "
                                                                                      + str(background_str)
                                                                                      )      
