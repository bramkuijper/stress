#!/usr/bin/env python3
import os, re, sys, math, datetime
import numpy as np
import socket as st
from scipy import optimize

# calculate optimal hormone levels in a non-autocorrelated
# environment where there are no lingering damage costs
class OptHorm:

    def __init__(self, ad, mu_bg, sP2NP, sNP2P, p_att, ap):
        self.ad = ad
        self.mu_bg = mu_bg
        self.pr_P = sNP2P / (sP2NP + sNP2P)
        self.p_att = p_att
        self.ap = ap

    def mort(self, h):
        return(self.mu_bg +\
                (1.0 - self.mu_bg) * self.pr_P * self.p_att * (1 - h**self.ap))

    def repr(self, h):
        return(1.0 - h**self.ad)

    def lrs(self, h):

        # reproductive success is (1 - mort) * repr
        # lifespan = 1/mort(h)
        return(-1*(1.0 - self.mort(h)) * self.repr(h) / self.mort(h))

    # print out hormone values and the corresponding lrs
    # function to check 
    def print_data(self):

        hval = list(np.linspace(0,1,50))

        the_str = "h;lrs;\n"

        for hval_i in hval:
            the_lrs = self.lrs(hval_i)
            the_str += f"{hval_i};{the_lrs}\n"

        print(the_str)

def logistic(x):
    return(1.0 / (1.0 + np.exp(-x)))

def logit(x):
    return(math.log(x/(1.0 - x)))


mu_feedback = [0.001]
mu_stress_cue_influx_combis = [[0.001,0.0]]
mu_influx = [0.001]
sdmu = [ 0.10]

# switch rate from P to NP
#s_P2NP = [[0.01,0.02],[ 0.04,0.08], [ 0.1,0.2]]
#s_NP2P = [[0.005,0.01], [ 0.02,0.04], [ 0.05, 0.1]]

s_P2NP = [ 0.1, 0.5, 0.9 ]
s_NP2P = [ 0.1, 0.5, 0.9 ]

#cue_P = [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ]
#cue_NP = [ 0.0, 0.2, 0.3, 0.4 ]

cue_P = [ 0.0 ]
cue_NP = [ 0.0 ]

# power of the damage cost function
ad = [ 2.0 ]

# power of the hormone level survival function
aP = [ 0.5, 0.9, 1.0, 1.5 ]

damage_decay = [ 1.0 ]
damage_due_to_hormone = [ 1.0 ]

# number of replicates
nrep = 3

ctr = 0

init_lfeedback = [-1.0]
init_lstress_influx = [-4.0]
init_lstress_baseline_influx = [-4.0]
init_lcue_influx = init_lstress_influx
init_linflux = [-1.5]

# attack probability 
p_att = [ 0.5 ]
pleiotropy = [ 0.0, 0.8 ]
maxtime = "200000"

# attack probability
mort_background = 0.001

exe = "./xstress"


host = st.gethostname()

background_jobs = False





if re.search("anthoxanthum",host) is None and background_jobs:
    background_jobs = False

background_str = " & " if background_jobs else ""

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
                                                        for init_lfeedback_i in init_lfeedback:
                                                            for init_lstress_influx_i in init_lstress_influx:
                                                                for init_lcue_influx_i in init_lcue_influx:
                                                                    for init_lstress_baseline_influx_i in init_lstress_baseline_influx:
                                                                        for pleiotropy_i in pleiotropy:

                                                                            ctr += 1
                                                                            print("echo " + str(ctr))

                                                                            base_name_i = base_name + "_"  + str(ctr)

                                                                            # optimal hormone object
                                                                            opt_h_obj = OptHorm(
                                                                                    ad=ad_i
                                                                                    ,mu_bg=mort_background
                                                                                    ,sP2NP=s_P2NP_i
                                                                                    ,sNP2P=s_NP2P_i
                                                                                    ,p_att=p_att_i
                                                                                    ,ap=aP_i)

#                                                                                opt_h_obj.print_data()

                                                                            sol = optimize.minimize_scalar(
                                                                                    fun=opt_h_obj.lrs
                                                                                    ,bounds=[1*10**(-5),1.0 - 1*10**(-5)]
                                                                                    ,method="Bounded"
                                                                                    )

                                                                            # optimal hormone level for 
                                                                            # this parameter combination
                                                                            h_opt = sol["x"]

                                                                            # now holding lfeedback value fixed
                                                                            # calculate value of influx
                                                                            influx = h_opt * logistic(init_lfeedback_i)

                                                                            # transform back to original logistic scale
                                                                            init_linflux_i = round(logit(influx),2)

                                                                            print(exe + " " 
                                                                                    + str(mu_feedback_i) + " " 
                                                                                    + str(mu_cue_influx_i) + " " 

                                                                                    + str(mu_stress_influx_i) + " " 
                                                                                    + str(mu_influx_i) + " " 

                                                                                    + str(sdmu_i) + " " 
                                                                                    + str(s_P2NP_i) + " " 

                                                                                    + str(s_NP2P_i) + " " 
                                                                                    + str(init_lfeedback_i) + " " 
                                                                                    + str(init_lcue_influx_i) + " " 
                                                                                    + str(init_lstress_influx_i) + " " 
                                                                                    + str(init_lstress_baseline_influx_i) + " " 
                                                                                    + str(init_linflux_i) + " " 
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
