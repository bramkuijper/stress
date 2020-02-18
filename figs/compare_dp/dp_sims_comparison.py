#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:45:19 2020

@author: bram
"""


import stress_panels
import pandas as pd
import numpy as np

def transform_hormone(row):
    return(row["mean_hormone"])/1000


dp_summary_file = "summary_dp.csv"
sims_summary_file = "summary_sims.csv"

translate_cols = {
    "t":"time"
    ,"pLeave":"sP2NP_1"
    ,"pArrive":"sNP2P_1"
    ,"pAttack":"p_att"
    ,"hormone":"mean_hormone"
    }

# read in the dp data summary and strip column values of white space
dp_data = pd.read_csv(filepath_or_buffer=dp_summary_file
                      ,sep="\t").rename(str.strip,axis="columns")

# translate column values
dp_data = dp_data.rename(
    columns=translate_cols
    )

dp_data["aP"] = 1.0
dp_data["ad"] = 1.5



dp_data["mean_hormone"] = dp_data.apply(transform_hormone,axis=1)

print(dp_data)

# read in simulation data
sim_data = pd.read_csv(filepath_or_buffer=sims_summary_file
                       ,sep=";")

# for testing purposes
params_panel_1 = {
    "aP": 1
    ,"ad":1.5
    ,"sP2NP_1":0.19
    ,"sNP2P_1":0.01}


stress_panels.stress_multipanel(
        param_array=np.array([[params_panel_1]])
        ,sim_data=sim_data
        ,dp_data=dp_data
        ,filename="test_plot.pdf"
        ,min_time_iter = 80
        ,max_time_iter = 175
        ,newzero = -100
        ,xlim=[-25,75]
        )