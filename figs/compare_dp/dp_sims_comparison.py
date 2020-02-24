#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:45:19 2020

@author: bram
"""


import stress_panels
import pandas as pd
import numpy as np



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
                      ,sep=";").rename(str.strip,axis="columns")

# translate column values
dp_data = dp_data.rename(
    columns=translate_cols
    )

dp_data["aP"] = 1.0
dp_data["ad"] = 1.5

print(dp_data.columns.values)

# read in simulation data
sim_data = pd.read_csv(filepath_or_buffer=sims_summary_file
                       ,sep=";")


params_panel_1 = {
    "aP": 1
    ,"ad":1.5
    ,"sP2NP_1":0.19
    ,"sNP2P_1":0.01}

params_panel_2 = {
    "aP": 1
    ,"ad":1.5
    ,"sP2NP_1":0.95
    ,"sNP2P_1":0.05}

stress_panels.stress_multipanel(
        param_array=np.array([[params_panel_1,params_panel_2]])
        ,title_array=[[r"Autocorrelation, $\rho = 0.8","No autocorrelation"]]
        ,sim_data=sim_data
        ,dp_data=dp_data
        ,filename="dp_vs_sims_main.svg"
        ,min_time_iter = 80
        ,max_time_iter =200 
        ,newzero = -100
        ,xlim=[-25,100]
        )

params_panel_1 = {
    "aP": 1
    ,"ad":1.5
    ,"sP2NP_1":0.09
    ,"sNP2P_1":0.01}

params_panel_2 = {
    "aP": 1
    ,"ad":1.5
    ,"sP2NP_1":0.9
    ,"sNP2P_1":0.1}

stress_panels.stress_multipanel(
        param_array=np.array([[params_panel_1,params_panel_2]])
        ,title_array=[[r"Autocorrelation, $\rho = 0.9","No autocorrelation"]]
        ,sim_data=sim_data
        ,dp_data=dp_data
        ,filename="dp_vs_sims_other.svg"
        ,min_time_iter = 80
        ,max_time_iter = 200
        ,newzero = -100
        ,xlim=[-25,100]
        )


params_panel_1 = {
    "aP": 1
    ,"ad":0.5
    ,"sP2NP_1":0.09
    ,"sNP2P_1":0.01}

params_panel_2 = {
    "aP": 1
    ,"ad":0.5
    ,"sP2NP_1":0.95
    ,"sNP2P_1":0.05}

stress_panels.stress_multipanel(
        param_array=np.array([[params_panel_1,params_panel_2]])
        ,title_array=[[r"Autocorrelation, $\rho = 0.9","No autocorrelation"]]
        ,sim_data=sim_data
        ,dp_data=dp_data
        ,filename="bang_bang.svg"
        ,min_time_iter = 80
        ,max_time_iter = 200
        ,newzero = -100
        ,xlim=[-25,100]
        )

