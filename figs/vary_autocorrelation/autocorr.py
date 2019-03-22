#!/usr/bin/env python3

# plot a phaseplot of pH and pD for different levels of dispersal

import pandas as pd
import functools
import itertools
import re, string, sys, os, math
import numpy as np
import subprocess
import matplotlib
matplotlib.use('pgf')
#matplotlib.use('Agg')


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

plt.style.use('base')

rcParams['axes.labelsize'] = 15

pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
         "\\usepackage{units}",         # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",  # unicode math setup
         r"\setmathfont{MyriadPro-Regular.otf}",
         r"\setmainfont{MyriadPro-Regular.otf}", # serif font via preamble
         ]
}

# add the latex stuff to the rcParams (parameters to plot the things)
matplotlib.rcParams.update(pgf_with_custom_preamble)


ctr = 0

# make stress iterations for this row of the data.frame
def stress_iterations(row, max_time=30, nrep=100):

    # get hormone level
    mean_hormone = row["mean_hormone"]
    var_hormone = row["var_hormone"]

    if var_hormone < 0:
        var_hormone = 0.0

    influx = row["mean_influx"]
    var_influx = row["var_influx"]

    if var_influx < 0:
        var_influx = 0.0

    stress_influx = row["mean_stress_influx"]
    var_stress_influx = row["var_stress_influx"]
    
    if var_stress_influx < 0:
        var_stress_influx = 0.0

    hormone_decay = row["mean_feedback"]
    var_hormone_decay = row["var_feedback"]
    
    if var_hormone_decay < 0:
        var_hormone_decay = 0.0

    zt_vals = np.random.normal(
            loc = mean_hormone
            ,scale = math.sqrt(var_hormone)
            ,size = nrep)

    decay_vals = np.random.normal(
            loc = hormone_decay
            ,scale = math.sqrt(var_hormone_decay)
            ,size = nrep
            )

    decay_vals = list(decay_vals)

    influx_vals = np.random.normal(
            loc = influx
            ,scale = math.sqrt(var_influx)
            ,size = nrep)

    influx_vals = list(influx_vals)

    stress_influx_vals = np.random.normal(
            loc = stress_influx
            ,scale = math.sqrt(var_stress_influx)
            ,size = nrep)

    stress_influx_vals = list(stress_influx_vals)

    iter_data = None

    # now go through all the sampled replicates
    for repl_i in range(nrep):

        # allocate a list of values of the stress level
        stress_t = [ 0 for t in range(max_time + 1) ]

        stress_t[0] = zt_vals[repl_i]

        stress_tplus1 = 0

        # iterate system for max_time timestpes
        for time_step in range(max_time):

            # get value for time t+1
            stress_t[time_step + 1] = stress_t[time_step] * (1.0 - decay_vals[repl_i]) + \
                    influx_vals[repl_i]

            if time_step == 10:
                stress_t[time_step + 1] = stress_t[time_step + 1] + stress_influx_vals[repl_i]

        current_iter_data = pd.DataFrame(
                    data=np.transpose(
                        np.array([list(range(max_time + 1))
                                    ,stress_t
                                    ,[ repl_i for i in range(max_time + 1) ]
                                ])
                        )
                    ,columns=["time","stress","replicate"]
                    )

        if iter_data is None:
            iter_data = current_iter_data
        else:
            iter_data = iter_data.append(current_iter_data, ignore_index=True)

    return(iter_data)

# the block function
def block(
        gs # gridspec
        ,row
        ,col
        ,data
        ,xlabel=None
        ,ylabel=None
        ,title=None
        ,plot_empty=False
        ,xtick_labels=False
        ,ytick_labels=False
        ,xlim=(-1,32)
        ,ylim=(0,200)
        ,ind_label=True):

    # get the counter for the letter indicator
    global ctr

    # initialize figure
    ax = plt.subplot(gs[row,col])
    
    # number of replicate simulations
    nrepl = data.shape[0]

    cmap = cm.get_cmap('tab10')

    for i in range(nrepl):

        # get data with a 40 timestep iteration
        # x 200 of stress iterations
        stress_iters = stress_iterations(data.iloc[i],max_time=30)

        replicate_id = list(stress_iters["replicate"].unique())

        for replicate_id_i in replicate_id:

            subset = stress_iters[stress_iters["replicate"] == replicate_id_i]

            ax.plot(subset["time"]
                    ,subset["stress"]
                    ,color=cmap.colors[i]
                    ,alpha=0.1
                    ,linewidth=0.5)

    # finalize the plot
    ax = finishblock(
            ax
            ,xlabel
            ,ylabel
            ,title
            ,plot_empty
            ,xtick_labels
            ,ytick_labels
            ,xlim=xlim
            ,ylim=ylim
            ,ind_label=ind_label)

    return(ax)


def finishblock(
        ax
        ,xlabel=None
        ,ylabel=None
        ,title=None
        ,plot_empty=False
        ,xtick_labels=False
        ,ytick_labels=False
        ,ylim=(-0.05,1.05)
        ,xlim=(-0.05,1.05)
        ,ind_label=True):

    global ctr

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    print(xlim)
    print(ylim)

    # do axis labeling
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if type(ylabel) == type("string"):
        ax.set_ylabel(ylabel=ylabel)
    
    if type(xlabel) == type("string"):
        ax.set_xlabel(xlabel=xlabel)

    if not xtick_labels:
        ax.xaxis.set_ticklabels([])

    if not ytick_labels:
        ax.yaxis.set_ticklabels([])
   
#    if ind_label:
#        ax.set_title(loc="left", 
#                label=string.ascii_uppercase[ctr], 
#                position=(0.0,1.02))

    ctr += 1
    
    if title is not None: 
        ax.set_title(
                label=title, 
                position=(0.9,0.1))


    return(ax)


# read in the data, where helping is conditional
data = pd.read_csv("../../data/summary_stress_autocorr.csv", sep=";")

cue_NP = []#list(data["cue_NP"].unique())
cue_P = [] #list(data["cue_P"].unique())

for cueNP_i in cue_NP:
    for cueP_i in cue_P:

        # initialize the figure
        fig = plt.figure(figsize=(30, 30))


        sp2np = [ 0, 0.1, 0.2, 0.3, 0.4 ] #list(data["sP2NP_2"].unique())
        snp2p = [ 0, 0.1, 0.2, 0.3, 0.4 ]#list(data["sNP2P_2"].unique())
        sp2np.sort()
        snp2p.sort()

        # set the graph's grid
        # with corresponding widths and heights
        widths = [ 1 for val in sp2np ]
        heights = [ 1 for val in snp2p ]


        # start gridspec object
        gs = gridspec.GridSpec(
                len(heights), 
                len(widths), 
                width_ratios=widths, 
                height_ratios=heights)

        for sp in sp2np:
            for snp in snp2p:

                data_subset = data[(data["sP2NP_2"] == sp)
                        & (data["sNP2P_2"] == snp) 
                        & (data["cue_P"] == cueP_i)
                        & (data["cue_NP"] == cueNP_i)
                        ]

                print(str(sp) + " " + str(snp))

                block(
                    gs
                    ,row=sp2np.index(sp)
                    ,col=snp2p.index(snp)
                    ,data = data_subset
                    ,xtick_labels = sp == max(sp2np)
                    ,ytick_labels = snp == max(snp2p)
                    ,title = r"$s_{\mathrm{p} \rightarrow \mathrm{np}} = " + str(sp) + "\n"
                            + r"$s_{\mathrm{np} \rightarrow \mathrm{p}} = " + str(snp)
                    )

        format = "pdf"
        graphname = "stress_curve_autocorrelation_cueNPi_" + str(cueNP_i) + "_cuePi" + str(cueP_i)
        graphname_pdf = graphname + ".pdf"
        graphname_svg = graphname + ".svg"
        plt.savefig(
                graphname_pdf 
                ,format="pdf"
                ,bbox_inches="tight"
        #        ,bbox_extra_artists=(x_label_txt,the_text,)
                ,transparent=True)

        if format == "svg":
            subprocess.call(["pdf2svg",graphname_pdf,graphname_svg])


        plt.close()


# initialize the figure
fig = plt.figure(figsize=(6, 6))
# set the graph's grid
# with corresponding widths and heights
widths = [ 1 ]
heights = [ 1 ]


# start gridspec object
gs = gridspec.GridSpec(
        len(heights), 
        len(widths), 
        width_ratios=widths, 
        height_ratios=heights)


sp = 0.1
snp = 0.1
cueP_i = 0.8
cueNP_i = 0.2

data_subset = data[(data["sP2NP_2"] == sp)
        & (data["sNP2P_2"] == snp) 
        & (data["cue_P"] == cueP_i)
        & (data["cue_NP"] == cueNP_i)
        ]


assert(data_subset.shape[0] > 0)

print(str(sp) + " " + str(snp))

block(
    gs
    ,row=0
    ,col=0
    ,data = data_subset
    ,xtick_labels =True
    ,ytick_labels = True
    ,title = r"$s_{\mathrm{p} \rightarrow \mathrm{np}} = " + str(sp) + "\n"
            + r"$s_{\mathrm{np} \rightarrow \mathrm{p}} = " + str(snp)
    )

format = "pdf"
graphname = "stress_curve_autocorrelation_cueNPi_" + str(cueNP_i) + "_cuePi" + str(cueP_i)
graphname_pdf = graphname + ".pdf"
graphname_svg = graphname + ".svg"
plt.savefig(
        graphname_pdf 
        ,format="pdf"
        ,bbox_inches="tight"
#        ,bbox_extra_artists=(x_label_txt,the_text,)
        ,transparent=True)
