#!/usr/bin/env python3



import stress_panels
import pandas as pd
import numpy as np

inv_u = True
inv_u_exp = False

dp_summary_file = "summary_dp.csv"
sims_summary_file = "summary_sims_standard.csv"
graph_name = "fourpanel_fig_base.pdf"
adval = 1.5
function_val = 1


if inv_u:
    dp_summary_file = "summary_dp_inverted_u.csv"
    sims_summary_file = "summary_sims_invu.csv"
#    sims_summary_file = "summary_sims_invu_long.csv"
    graph_name = "fourpanel_fig_invu.pdf"
    adval = 2


    if inv_u_exp:
        dp_summary_file = "summary_dp_inverted_u_exp.csv"
        sims_summary_file = "summary_sims_invu_long.csv"
        graph_name = "fourpanel_fig_invu_exp.pdf"
        function_val = 2

translate_cols = {
    "t":"time"
    ,"pLeave":"sP2NP_1"
    ,"pArrive":"sNP2P_1"
    ,"pAttack":"p_att"
    ,"hormone":"mean_hormone"
    }

dp_data = pd.read_csv(filepath_or_buffer=dp_summary_file
                      ,sep=";").rename(str.strip,axis="columns")

# translate column values
dp_data = dp_data.rename(
    columns=translate_cols
    )

dp_data["aP"] = 1.0
dp_data["ad"] = adval


# read in simulation data
sim_data = pd.read_csv(filepath_or_buffer=sims_summary_file
                       ,sep=";")

sim_data["ad"] = adval

if inv_u:
    sim_data = sim_data.loc[sim_data["fecundity_function"] == function_val]

params_panel_1 = {
    "aP": 1
    ,"ad":adval
    ,"sP2NP_1":0.95
    ,"sNP2P_1":0.05}

params_panel_2 = {
    "aP": 1
    ,"ad":adval
    ,"sP2NP_1":0.665
    ,"sNP2P_1":0.035}

params_panel_3 = {
    "aP": 1
    ,"ad":adval
    ,"sP2NP_1":0.9
    ,"sNP2P_1":0.1}

params_panel_4 = {
    "aP": 1
    ,"ad":adval
    ,"sP2NP_1":0.095
    ,"sNP2P_1":0.005}

in_panel_text1 = dict(x=0.7
        ,y=0.1
        ,label="Risk: 0.05" + "\n" + "Autocorrelation: 0")

in_panel_text2 = dict(x=0.5
        ,y=0.5
        ,label="Risk: 0.05" + "\n" + "Autocorrelation: 0.3")

in_panel_text3 = dict(x=0.7
        ,y=0.1
        ,label="Risk: 0.1" + "\n" + "Autocorrelation: 0")

in_panel_text4 = dict(x=0.5
        ,y=0.5
        ,label="Risk: 0.1" + "\n" + "Autocorrelation: 0.9")


text_array = np.array([[in_panel_text1, in_panel_text2],[in_panel_text3,in_panel_text4]])

stress_panels.stress_multipanel(
        param_array=np.array([[params_panel_1,params_panel_2],[params_panel_3, params_panel_4]])
        ,title_array=[[r"No autocorrelation","Autocorrelation"],["",""]]
        ,label_array=text_array
        ,sim_data=sim_data
        ,dp_data=dp_data
        ,filename=graph_name
        ,min_time_iter = 80
        ,max_time_iter =200 
        ,newzero = -100
        ,xlim=[-25,100]
        )


