#!/usr/bin/env python3
from pythontools import multipanel
import itertools
import argparse
import pandas as pd
import numpy as np
import sys, re, os
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.use("pgf")

# font hassle 
# if one gets errors with font not found
# just comment this custom preamble stuff out
# it is just that the default fonts in matplotlib
# suck a bit
pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
        r"\usepackage{units}",         # load additional packages
        r"\usepackage{mathspec}",         # load additional packages
        r"\setmainfont[" +\
            "UprightFont = * ," +\
            "ItalicFont = *Oblique," +\
            "BoldFont = *Bold ," +\
            "Extension = .otf]{FreeSans}",
        r"\setmathsfont(Digits,Latin,Greek)[" +\
            "UprightFont = * ," +\
            "ItalicFont = *Oblique," +\
            "BoldFont = *Bold ," +\
            "Extension = .otf]{FreeSans}",
        r"\setmathrm[" +\
            "UprightFont = * ," +\
            "ItalicFont = *Oblique," +\
            "BoldFont = *Bold ," +\
            "Extension = .otf]{FreeSans}",
         ]
}

mpl.rcParams.update(pgf_with_custom_preamble)

# some further parameter specifications
mpl.rcParams["axes.labelsize"] = 18
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["legend.fontsize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"


# set up command-line argument parsing
# to specify whether we want juvenile phenotypes
# plotted or not (this makes it easier to call the 
# same script multiple times
parser = argparse.ArgumentParser()

# specify output file name
parser.add_argument('-o', default="example_simulation.pdf")
parser.add_argument('-i'
        ,nargs="+" # provide multiple simulation files
        ,help="file containing single simulation over time")
args = vars(parser.parse_args())


def skip_parameters_front(filename):

    """
    opens a simulation file, reads outptu and looks for the line where
    the data starts and the parameter listing ends

    Parameters:
    -----------
    filename: str
        the name fo the file containing a single simulation

    Returns:
    int
        the line number where the data starts

    """

    f = open(filename)
    fl = f.readlines()
    f.close()
    linerange = list(range(0,len(fl)))

    for line_i in fl:
        if re.search("generation",line_i) is not None:
            return(line_i)



def get_parameter_line(filename):

    """
    opens a simulation file, reads output and looks for the line
    where the data ends and the parameter listing begins

    Parameters:
    -----------
    filename: str
        the name of the file containing a single simulation

    Returns: 
    int
        the line number where the data ends
    """

    # open file read lines and close
    f = open(filename)
    fl = f.readlines()
    f.close()

    # get a reversed list of line numbers
    # (we are searching from the bottom of the file
    # upwards, saving us having to go through the whole of the file)
    linerange = list(range(0,len(fl)))
    linerange.reverse()

    for line_i in linerange:

        # look for the first line which starts with a number
        # this is a line containing simulation data and hence
        # is not the parameter listing
        if re.search("^\d",fl[line_i]) is not None:
            return(line_i)


##### get the data from the simulation and select subsets #####

# list of simulation files
simulation_folder = args["i"]

datasets = []

file_list = os.listdir(path=simulation_folder)

for file_i in file_list:

    if re.search("sim_.*",file_i) is not None:
        startline = skip_parameters_front(file_i)
        set_i = pd.read_csv(file_i,sep=";",skiprows=parline-1)

        datasets += [set_i]


data_generation_interval = datasets[0]["generation"].max() / 6


#### initialize figure ####

#panels
# 1.stress influx and baseline influx
# 2. feedback 
# 3. elevation
# 4. environment


# give the figure twice as many columns as omega values
# as each block contains one panel for Sonora and one for Catalina
the_fig = multipanel.MultiPanel(
        panel_widths=[1,0.1,1]
        ,panel_heights=[1,1,1,1]
        ,filename=args["o"] # value of the output file
        ,width=15
        ,height=15
        )


#### phenotype ####
the_axis = the_fig.start_block(
            row=0
            ,col=0)

# dictionary with colors labels etc
phen_dict = {
        "juv": {"color":"#ff0000", "label":r"$\bar{z}_{\mathrm{lv}}$"}
        ,"ad": {"color":"#0000ff", "label":r"$\bar{z}_{\mathrm{ad}}$"}
}

stress_influx_color = "#ff0069"
influx_color = "#00bdff"

for i, dataset_i in enumerate(datasets):

    the_axis.plot(
            dataset_i["generation"]
            ,dataset_i["mean_stress_influx"]
            ,color=stress_influx_color
            ,label="_nolabel" if i > 0 else "Stress influx"
            ,linewidth=1
            )

    the_axis.plot(
            dataset_i["generation"]
            ,dataset_i["mean_influx"]
            ,color=influx_color
            ,label="_nolabel" if i > 0 else "Baseline influx"
            ,linewidth=1
            )

the_axis.legend()

ylim = [-0.05,1.05]


# end the figure
# see multipanel module
the_fig.end_block(
        the_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.25
        ,xticks=False
        ,yticks=True
        ,title=""
        ,ylabel=r"Influx"
        ,xlabel=""
        ,loc_title=True
        )

#### genetic elevation ####
feedback_axis = the_fig.start_block(
            row=1
            ,col=0)

feedback_color = "#53128d"

for i, dataset_i in enumerate(datasets):
    
    feedback_axis.plot(
            dataset_i["generation"]
            ,1.0 - dataset_i["mean_feedback"]
            ,color=feedback_color
            ,label="_nolabel" if i > 0 else "Feedback"
            ,linewidth=1
            )

ylim = [-0.05,1.05]
# end the figure
# see multipanel module
the_fig.end_block(
        feedback_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.25
        ,xticks=True
        ,yticks=True
        ,title=""
        ,ylabel=r"Feedback, $a$"
        ,loc_title=True
        )

#### hormone ####
hormone_axis = the_fig.start_block(
            row=0
            ,col=2)

hormone_color = "#5c8d12"

for i, dataset_i in enumerate(datasets):

    hormone_axis.plot(
            dataset_i["generation"]
            ,dataset_i["mean_hormone"] 
            ,color=hormone_color
            ,alpha=1.0
            ,linewidth=1)

ylim = [ -0.05, 1.05]

# end the figure
# see multipanel module
the_fig.end_block(
        hormone_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.25
        ,xticks=False
        ,yticks=True
        ,title=""
        ,ylabel=r"Hormone, $z$"
        ,xlabel=""
        ,loc_title=True
        )


#### first within-generational plasticity ####
damage_axis = the_fig.start_block(
            row=0
            ,col=3)

damage_color = "#8d4112"

for i, dataset_i in enumerate(datasets):

    damage_axis.plot(
            data["generation"]
            ,data["mean_damage"] 
            ,color=damage_color
            ,alpha=1.0
            ,linewidth=1)

ylim = [ -0.05, 1.05 ] 

the_fig.end_block(
        mat_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.25
        ,xticks=True
        ,yticks=True
        ,title=""
        ,ylabel=r"Damage"
        ,xlabel=r"Time, $t$"
        ,loc_title=True
        )

the_fig.close(tight=True)
