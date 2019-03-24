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


#########################################
#           read in the data
#########################################

line_num_params = find_parameter_linenum(sys.argv[1])

# get the data
dat = pd.read_csv(sys.argv[1],
        nrows=line_num_params-1,
        sep=';',
        index_col=False)

#########################################
#           make figure
#########################################

# initialize the figure
fig = plt.figure(figsize=(10,15))

# generate the grid of the graph
# see: 
widths = [ 1 ]
heights = [ 1, 1, 1, 1, 1, 1, 1]
numrows = len(heights)
numcols  = len(widths)

# make the grid
gs = gridspec.GridSpec(
        nrows=numrows,
        ncols=numcols,
        width_ratios=widths,
        height_ratios=heights)


# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[0,0])

ax.plot(
        dat["generation"],
        dat["mean_feedback"],
        label="feedback",
        color="red"
        )

ax.fill_between(
        x=dat["generation"],
        y1=dat["mean_feedback"] + dat["var_feedback"].apply(np.sqrt),
        y2=dat["mean_feedback"] - dat["var_feedback"].apply(np.sqrt),
        label="_nolabel",
        color="red",
        alpha=0.1,
        linewidth=0
        )

ax.set_ylabel("Mean" + "\n" +  r"Feedback$\pm$SD")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)


# plot influxes 
ax = plt.subplot(gs[1,0])

ax.set_ylabel(r"Influx")

ax.plot(
        dat["generation"],
        dat["mean_influx"],
        label="Influx",
        color="blue")

ax.fill_between(
        x=dat["generation"],
        y1=dat["mean_influx"] + dat["var_influx"].apply(np.sqrt),
        y2=dat["mean_influx"] - dat["var_influx"].apply(np.sqrt),
        label="_nolabel",
        color="blue",
        alpha=0.1,
        linewidth=0
        )

ax.plot(
        dat["generation"],
        dat["mean_stress_influx"],
        label="Stress influx",
        color="red")

ax.fill_between(
        x=dat["generation"],
        y1=dat["mean_stress_influx"] + dat["var_stress_influx"].apply(np.sqrt),
        y2=dat["mean_stress_influx"] - dat["var_stress_influx"].apply(np.sqrt),
        label="_nolabel",
        color="red",
        alpha=0.1,
        linewidth=0
        )

ax.legend()

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

####### plot damage #######
ax = plt.subplot(gs[2,0])

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
ax = plt.subplot(gs[3,0])

ax.plot(
        dat["generation"],
        dat["mean_hormone"],
        label="hormone")

ax.set_ylabel(r"Hormone")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

####### plot variances #######
ax = plt.subplot(gs[4,0])

ax.plot(
        dat["generation"],
        dat["var_feedback"],
        label="Feedback")

ax.plot(
        dat["generation"],
        dat["var_stress_influx"],
        label="Stress influx")

ax.plot(
        dat["generation"],
        dat["var_influx"],
        label="Influx")

ax.plot(
        dat["generation"],
        dat["var_hormone"],
        label="Hormone")

ax.plot(
        dat["generation"],
        dat["var_damage"],
        label="Damage")

ax.legend()
ax.set_ylabel(r"Var")

ax.set_xlabel(r"Generation")

####### plot stress response curve over time #######

ax = plt.subplot(gs[5,0])

zt = float(dat["mean_hormone"][-1:])
var_zt = float(dat["var_hormone"][-1:])

feedback = float(dat["mean_feedback"][-1:])
var_feedback = float(dat["var_feedback"][-1:])

influx = float(dat["mean_influx"][-1:])
var_influx = float(dat["var_influx"][-1:])

stress_influx = float(dat["mean_stress_influx"][-1:])
var_stress_influx = float(dat["var_stress_influx"][-1:])

nrep = 2000

print(zt)

for i in range(0,nrep):

    zt_val = np.random.normal(
            loc=zt, 
            scale=var_zt,
            size=1)

    feedback_val = np.random.normal(
            loc=feedback, 
            scale=np.sqrt(var_feedback),
            size=1)

    influx_val = np.random.normal(
            loc=influx, 
            scale=var_influx,
            size=1)

    stress_influx_val = np.random.normal(
            loc=stress_influx, 
            scale=var_stress_influx,
            size=1)

    time = list(range(0,30))

    stress = [ 0 for x in time ]
    stress[0] = zt_val

    for t in time[0:-1]:
        stress[t+1] = stress[t] * (1 - feedback_val) + influx_val

        if t == 10:
            stress[t+1] += stress_influx_val

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

ax.plot(
        time[5:]
        ,stress[5:]
        ,label="_nolabel"
        ,linewidth=2
        ,color="black")


ax.set_ylim(0, zt*2)

ax.set_ylabel(r"Stress response" + "\n" + r"to stimulus at $t=10$")
ax.set_xlabel(r"Time")


ax = plt.subplot(gs[6,0])

ax.plot(
        dat["generation"],
        dat["var_feedback"],
        label="Feedback")


format = "pdf"

filename = os.path.join(
        os.path.dirname(sys.argv[1]),
        "graph_" + os.path.basename(sys.argv[1]) + "." + format
        )

print(type(filename))

plt.savefig(
        fname=filename,
        format=format, 
        bbox_inches="tight")
