#!/usr/bin/env python3

# plot a single individual based simulation
# by doing ./plot_simulation.py $filename

# load libraries 
# (mainly pandas, numpy and matplotlib)
import pandas as pd
import itertools
import subprocess 
import math
import argparse
import numpy as np
import sys, re, os.path
import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from pathlib import PurePath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

# some stuff to render fonts in graphs
rcParams['axes.labelsize'] = 15
rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# some stuff to render fonts in graphs 
# using latex
# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

#########################################
#           check where data ends
#########################################
def find_parameter_linenum(filename):

    # get all the lines
    f = open(filename)
    fl = f.readlines()
    f.close()

    # find the first line from below that
    # starts with a number. That is the last
    # line of the data
    frange = list(range(0,len(fl)))
    frange.reverse()

    # ok loop over all lines (backwards)
    for line_num in frange:
        if re.match("^\d.*",fl[line_num]) is not None:
            return(line_num+1)

    return(len(fl))


def check_params_first(file, regex=r"^generation"):

    """
    gets first line of dataset if parameters
    are printed before the actual data

    Parameters:
    -----------
    file:   str
        name of the data file

    Returns:
    --------
    linenum: int
        line number in file where data header starts
    """

    f = open(file)
    fl = f.readlines()
    f.close()

    for linenum, line in enumerate(fl):

        if re.search(regex, line) is not None:

            return(linenum)

#########################################
#           read in the data
#########################################
params_at_start = check_params_first(sys.argv[1])

assert(params_at_start is not None)

line_num_params = find_parameter_linenum(sys.argv[1])

print("number of rows to skip: "  + str(params_at_start))
print("endline "  + str(line_num_params))

# get the data
dat = pd.read_csv(sys.argv[1],
    nrows=line_num_params-1-params_at_start,
    sep=';',
    skiprows=params_at_start,
    index_col=False)

#########################################
#           make figure
#########################################

# initialize the figure
fig = plt.figure(figsize=(10,18))

# generate the grid of the graph
# see: 
widths = [ 1 ]
heights = [ 1 for i in range(0,9) ]
numrows = len(heights)
numcols  = len(widths)

# make the grid
gs = gridspec.GridSpec(
    nrows=numrows,
    ncols=numcols,
    width_ratios=widths,
    height_ratios=heights)

row_ctr = 0 

ax = plt.subplot(gs[row_ctr,0])

ax.plot(
    dat["generation"],
    dat["mean_feedback"],
    label="feedback",
    color="red"
    )

ax.fill_between(
    x=dat["generation"]
    ,y1=dat["mean_feedback"] + dat["sd_feedback"]
    ,y2=dat["mean_feedback"] - dat["sd_feedback"]
    ,label="_nolabel"
    ,color="red"
    ,alpha=0.1
    ,linewidth=0
    )

ax.set_ylabel(r"Feedback, $f$" + "\n" + r"where $h_{t+1} = (1-f) h_{t}$")

ax.tick_params(
    axis="x",
    which="both",
    labelbottom=False)


row_ctr += 1


# plot influxes 
ax = plt.subplot(gs[row_ctr,0])

ax.set_ylabel(r"Influx")

ax.fill_between(
    x=dat["generation"],
    y1=dat["mean_influx"] + dat["sd_influx"],
    y2=dat["mean_influx"] - dat["sd_influx"],
    label="_nolabel",
    color="blue",
    alpha=0.1,
    linewidth=0
    )

ax.plot(
    dat["generation"],
    dat["mean_influx"],
    label="Influx",
    color="blue")


ax.fill_between(
    x=dat["generation"],
    y1=dat["mean_stress_influx"] + dat["sd_stress_influx"],
    y2=dat["mean_stress_influx"] - dat["sd_stress_influx"],
    label="_nolabel",
    color="red",
    alpha=0.1,
    linewidth=0
    )

ax.plot(
    dat["generation"],
    dat["mean_stress_influx"],
    label="Stress influx",
    color="red")

if "mean_cue_influx" in dat.columns.values:

    ax.fill_between(
        x=dat["generation"],
        y1=dat["mean_cue_influx"] + dat["sd_cue_influx"],
        y2=dat["mean_cue_influx"] - dat["sd_cue_influx"],
        label="_nolabel",
        color="green",
        alpha=0.1,
        linewidth=0
        )

    ax.plot(
        dat["generation"],
        dat["mean_cue_influx"],
        label="Cue influx",
        color="green")

ax.legend()

ax.tick_params(
    axis="x",
    which="both",
    labelbottom=False)

####### plot damage #######
row_ctr += 1
ax = plt.subplot(gs[row_ctr,0])

ax.plot(
    dat["generation"],
    dat["mean_damage"],
    label="damage")

ax.set_ylabel(r"Damage")

ax.tick_params(
    axis="x",
    which="both",
    labelbottom=False)


####### plot hormones #######
row_ctr += 1 
ax = plt.subplot(gs[row_ctr,0])

ax.plot(
    dat["generation"],
    dat["mean_hormone"],
    label="hormone")

ax.set_ylabel(r"Hormone")

ax.tick_params(
    axis="x",
    which="both",
    labelbottom=False)

####### plot stdevs #######
row_ctr += 1
ax = plt.subplot(gs[row_ctr,0])

ax.plot(
    dat["generation"],
    dat["sd_feedback"],
    label="Feedback")

ax.plot(
    dat["generation"],
    dat["sd_stress_influx"],
    label="Stress influx")

ax.plot(
    dat["generation"],
    dat["sd_influx"],
    label="Influx")

ax.plot(
    dat["generation"],
    dat["sd_hormone"],
    label="Hormone")

ax.plot(
    dat["generation"],
    dat["sd_damage"],
    label="Damage")

ax.legend()
ax.set_ylabel(r"SD")

ax.set_xlabel(r"Generation")

####### plot stress response curve over time #######

row_ctr += 1
ax = plt.subplot(gs[row_ctr,0])

zt = float(dat["mean_hormone"][-1:])
sd_zt = float(dat["sd_hormone"][-1:])

feedback = float(dat["mean_feedback"][-1:])
sd_feedback = float(dat["sd_feedback"][-1:])

influx = float(dat["mean_influx"][-1:])
sd_influx = float(dat["sd_influx"][-1:])

stress_influx = float(dat["mean_stress_influx"][-1:])
sd_stress_influx = float(dat["sd_stress_influx"][-1:])

nrep = 1000

for i in range(0,nrep):

    zt_val = np.random.normal(
            loc=zt, 
            scale=sd_zt,
            size=1)

    feedback_val = np.random.normal(
            loc=feedback, 
            scale=sd_feedback,
            size=1)

    influx_val = np.random.normal(
            loc=influx, 
            scale=sd_influx,
            size=1)

    stress_influx_val = np.random.normal(
            loc=stress_influx, 
            scale=sd_stress_influx,
            size=1)

    time = list(range(0,30))

    stress = [ 0 for x in time ]
    stress[0] = zt_val

    for t in time[0:-1]:
        stress[t+1] = stress[t] * (1 - feedback_val) + influx_val

        if t == 10:
            stress[t+1] += stress_influx_val

        if stress[t+1] > 1.0:
            stress[t+1] = 1.0

    # plot each individual stress curve
    ax.plot(
            time[5:]
            ,stress[5:]
            ,label="_nolabel"
            ,alpha=0.04
            ,color="blue"
            )

# now calculate the mean stress level
time = list(range(0,30))

stress = [ 0 for x in time ]
stress[0] = zt

for t in time[0:-1]:
    stress[t+1] = stress[t] * (1 - feedback) + influx

    if t == 10:
        stress[t+1] += stress_influx
        
    if stress[t+1] > 1.0:
        stress[t+1] = 1.0

ax.plot(
        time[5:]
        ,stress[5:]
        ,label="_nolabel"
        ,linewidth=2
        ,color="black")


ax.set_ylim([-0.05, 1.05])

ax.set_ylabel(r"Stress response" + "\n" + r"to stimulus at $t=10$")
ax.set_xlabel(r"Time")


#### mortality  ####
row_ctr += 1
ax = plt.subplot(gs[row_ctr,0])

ax.plot(dat["generation"]
        ,dat["prop_dead_P"]
        ,label="P"
        )

ax.plot(dat["generation"]
        ,dat["prop_dead_NP"]
        ,label="NP"
        )

ax.set_ylim([-0.05,1.05])
ax.legend()

ax.set_ylabel(r"Prob mort")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)
        
#### frequency of P individuals ####
ax = plt.subplot(gs[7,0])

ax.plot(dat["generation"]
        ,dat["freq_P"])

ax.set_ylim([-0.05,1.05])
ax.set_ylabel(r"Freq(P)")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)


#### fecundity ####

ax = plt.subplot(gs[8,0])

ax.plot(dat["generation"]
        ,dat["mean_fecundity_P"]
        ,label="P"
        )

ax.plot(dat["generation"]
        ,dat["mean_fecundity_NP"]
        ,label="NP"
        )
ax.legend()

ax.set_ylim([-0.05,1.05])
ax.set_ylabel(r"Fecundity")


format = "png"

filename = os.path.join(
        os.path.dirname(sys.argv[1]),
        "graph_" + os.path.basename(sys.argv[1]) + "." + format
        )

print("finalizing plot")
plt.savefig(
        fname=filename,
        format=format, 
        bbox_inches="tight")
