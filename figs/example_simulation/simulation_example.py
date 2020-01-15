#!/usr/bin/env python3
from pythontools import multipanel
import itertools
import argparse
import pandas as pd
import numpy as np
import sys, re
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
simulation_files = args["i"]
line_number_end = get_parameter_line(simulation_file)

data = pd.read_csv(
        simulation_file
        ,sep=";"
        ,nrows=line_number_end)

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
phen_axis = the_fig.start_block(
            row=0
            ,col=0)
    
# dictionary with colors labels etc
phen_dict = {
        "juv": {"color":"#ff0000", "label":r"$\bar{z}_{\mathrm{lv}}$"}
        ,"ad": {"color":"#0000ff", "label":r"$\bar{z}_{\mathrm{ad}}$"}
}

for key, value in phen_dict.items():

    mean_name = "meanz_" + key
#    var_name = "varz_" + key

    color = value["color"]
    label = value["label"]

#    phen_axis.fill_between(
#            x=data["generation"]
#            ,y1=data[mean_name] + np.sqrt(data[var_name])
#            ,y2=data[mean_name] - np.sqrt(data[var_name])
#            ,color=color
#            ,alpha=0.5
#            ,label="_nolabel"
#            ,linewidth=0.5)

    phen_axis.plot(
            data["generation"]
            ,data[mean_name] 
            ,color=color
            ,label=label
            ,alpha=1.0
            ,linewidth=1)


phen_axis.legend()

ylim = [ -2, 2]


# end the figure
# see multipanel module
the_fig.end_block(
        phen_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 1
        ,xticks=False
        ,yticks=True
        ,title=""
        ,ylabel=r"Phenotypes, $z$"
        ,xlabel=""
        ,loc_title=True
        )

#### genetic elevation ####
genetic_axis = the_fig.start_block(
            row=1
            ,col=0)

mean_name = "meang"
var_name = "varg"

color = "#008c22"

genetic_axis.fill_between(
        x=data["generation"]
        ,y1=data[mean_name] + np.sqrt(data[var_name])
        ,y2=data[mean_name] - np.sqrt(data[var_name])
        ,color=color
        ,alpha=0.5
        ,label="_nolabel"
        ,linewidth=0.5)

genetic_axis.plot(
        data["generation"]
        ,data[mean_name] 
        ,color=color
        ,label=label
        ,alpha=1.0
        ,linewidth=1)

ylim = [-0.25,0.25]
# end the figure
# see multipanel module
the_fig.end_block(
        genetic_axis 
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.1
        ,xticks=True
        ,yticks=True
        ,title=""
        ,ylabel=r"Elevation, $a$"
        ,xlabel=r"Generation, $t$"
        ,loc_title=True
        )
#### within-generational plasticity ####
plasticity_axis = the_fig.start_block(
            row=0
            ,col=2)

# dictionary with colors labels etc
plast_dict = {
        "juv": {"color":"#0bd9ff", "label":r"$\bar{b}_{\mathrm{lv}}$"}
        ,"ad": {"color":"#ff0068", "label":r"$\bar{b}_{\mathrm{ad}}$"}
}

for key, value in plast_dict.items():

    mean_name = "meanb_" + key
    var_name = "varb_" + key

    color = value["color"]
    label = value["label"]

    plasticity_axis.fill_between(
            x=data["generation"]
            ,y1=data[mean_name] + np.sqrt(data[var_name])
            ,y2=data[mean_name] - np.sqrt(data[var_name])
            ,color=color
            ,alpha=0.5
            ,label="_nolabel"
            ,linewidth=0.5)

    plasticity_axis.plot(
            data["generation"]
            ,data[mean_name] 
            ,color=color
            ,label=label
            ,alpha=1.0
            ,linewidth=1)


plasticity_axis.legend()

ylim = [ 0, 1]


# end the figure
# see multipanel module
the_fig.end_block(
        plasticity_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.25
        ,xticks=False
        ,yticks=True
        ,title=""
        ,ylabel=r"Plasticity, $b$"
        ,xlabel=""
        ,loc_title=True
        )


#### first within-generational plasticity ####
mat_axis = the_fig.start_block(
            row=1
            ,col=2)

# dictionary with colors labels etc
mat_dict = {
        "juv": {"color":"#0bd9ff", "label":r"$\bar{m}_{\mathrm{lv}}$"}
        ,"ad": {"color":"#ff0068", "label":r"$\bar{m}_{\mathrm{ad}}$"}
}

for key, value in mat_dict.items():

    mean_name = "meanm_m_" + key
    var_name = "varm_m_" + key

    color = value["color"]
    label = value["label"]

    mat_axis.fill_between(
            x=data["generation"]
            ,y1=data[mean_name] - np.sqrt(data[var_name])
            ,y2=data[mean_name] + np.sqrt(data[var_name])
            ,color=color
            ,alpha=0.5
            ,label="_nolabel"
            ,linewidth=0.5)

    mat_axis.plot(
            data["generation"]
            ,data[mean_name] 
            ,color=color
            ,label=label
            ,alpha=1.0
            ,linewidth=1)


mat_axis.legend()

ylim = [ -0.25, 1]

the_fig.end_block(
        mat_axis
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 0.25
        ,xticks=False
        ,yticks=True
        ,title=""
        ,ylabel=r"Maternal effect, $m$"
        ,xlabel=""
        ,loc_title=True
        )



#### environment ####
envt_axis = the_fig.start_block(
            row=2
            ,col=2)

color = "#ac7c00"

envt_axis.plot(
        data["generation"]
        ,data["epsilon"] 
        ,color=color
        ,alpha=1.0
        ,linewidth=0.5)

ylim = [ -2.5,2.5]

# end the figure
# see multipanel module
the_fig.end_block(
        envt_axis 
        ,ylim=ylim
        ,y_ticks_minor = 5
        ,x_ticks_minor = 5
        ,x_ticks_major_multiple = data_generation_interval
        ,y_ticks_major_multiple = 1
        ,xticks=True
        ,yticks=True
        ,title=""
        ,ylabel=r"Standardized temperature, $\varepsilon$"
        ,xlabel=r"Generation, $t$"
        ,loc_title=True
        )

the_fig.close(tight=True)
