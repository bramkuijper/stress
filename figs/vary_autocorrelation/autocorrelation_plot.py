#!/usr/bin/env python3

# plots the first generation's outcome
from pythontools import multipanel
import pandas as pd
from matplotlib import cm
import os.path
import sys

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

mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["axes.labelpad"] = 16
mpl.rcParams["svg.fonttype"] = "none"



### read in the data
the_data = pd.read_csv("summary_vary_autocorrelation_new.csv",sep=";")

#the_data = the_data.loc[(the_data["sP2NP_1"].isin([0.01,0.1,0.2,0.5]))
#        & (the_data["sNP2P_1"].isin([0.01,0.1,0.2,0.5]))]
#
assert(the_data.shape[0] > 0)

# TODO select only one particular row or so 


iterfolder = "hpcbatch_14_01_2020_225417_iter/"
############## auxiliary functions ##############
print("it's go...")

# simple function to use a dict to query a dataframe
def make_query(the_dict):

    qry_list = []

    for key, value in the_dict.items():
        qry_list += [str(key) + " == '" + str(value) + "'"]

    qry = " & ".join(qry_list)

    return(qry)

# from a summary data frame obtain filenames
# and extract iteration data from those
def get_iter_data(df, iter_dir):

    # get list of filenames
    file_list = list(df["file"].unique())

    iter_df = None

    # now open files and extract data
    # for each combination of parameters there are likely
    # multiple 
    for replicate, file_name in enumerate(file_list):

        # get the filename of the filename containing
        # the iterations
        file_base_name = os.path.basename(file_name)

        iter_file_name = os.path.join(iter_dir, file_base_name + "iters")

        try:
            subdf = pd.read_csv(iter_file_name,sep=";")

            subdf["replicate"] = replicate

            if iter_df is None:
                iter_df = subdf
            else:
                iter_df = iter_df.append(subdf, ignore_index=True)

        except:
            print("meh file " + iter_file_name + " won't exist")

    return(iter_df)

######## general plotting function ########

def stress_multipanel(
        parameter_dict
        ,filename
        ,max_time_iter = 150
        ):
   
    # use the dictionary to make a subset of the 
    # data
    query_string = make_query(parameter_dict)

    print(query_string)
    subset = the_data.query(query_string)

    assert(subset.shape[0] > 0)

    # get switching rates, which we use for the
    # individual panels
    p_leave_u = list(subset["sP2NP_1"].unique())
    p_arrive_u = list(subset["sNP2P_1"].unique())

    # initialize multipanel fig
    the_fig = multipanel.MultiPanel(
            panel_widths=[1 for x in p_leave_u]
            ,panel_heights=[1 for y in p_arrive_u]
            ,filename=filename
            ,width=5 * len(p_leave_u)
            ,height=5 * len(p_arrive_u)
            ,hspace = 0.3
            )

    # sort
    p_leave_u.sort()

    # we want the reverse sorting done
    # for the y axis, which should descend
    # from high to low
    p_arrive_u.sort()
    p_arrive_u.reverse()

    # get colors
    the_color_map = cm.get_cmap("tab10")

    # now loop over rows and columns
    for row_i, p_arrive_i in enumerate(p_arrive_u):
        for col_i, p_leave_i in enumerate(p_leave_u):

            # subset selection
            subset_subset = subset[
                    (subset["sP2NP_1"] == p_leave_i)
                    & (subset["sNP2P_1"] == p_arrive_i)
                    ]

            # obtain the iteration file data
            iteration_data = get_iter_data(subset_subset,iterfolder)

            ######## make each panel ########
            the_axis = the_fig.start_block(
                row=row_i
                ,col=col_i)

            # check whether there is iteration data for this
            # data set

            # if yes, plot it
            if iteration_data is not None:

                iteration_data = iteration_data[
                        iteration_data["time"] < max_time_iter 
                            ]


                # get list of unique replicates per parameter combination
                replicates = list(iteration_data["replicate"].unique())
                # get list of unique individuals per replicate
                individuals = list(iteration_data["individual"].unique())[0:10]

                print("number replicates per plot: " + str(len(replicates)))

                for replicate_i in replicates:
                    for individual_i in individuals:

                        # select the subset for individual i and replicate j
                        iteration_subset = iteration_data[
                                (iteration_data["replicate"] == replicate_i)
                                &  (iteration_data["individual"] == individual_i)]

                        the_axis.plot(iteration_subset["time"]
                                ,iteration_subset["stress"]
                                ,color=the_color_map.colors[replicate_i]
                                ,linewidth=0.5
                                ,alpha=0.5
                                ,label="_nolabel")

            ylab = ""

            # end the figure
            the_fig.end_block(
                    the_axis
                    ,ylim=[0,1]
                    ,y_ticks_minor = 1
                    ,x_ticks_minor = 1
                    ,x_ticks_major_multiple = 25
                    ,y_ticks_major_multiple = 0.25
                    ,xticks=True
                    ,yticks=col_i == 0
                    ,title=r"$\lambda_{P\rightarrow NP} = " + str(p_leave_i) + r"$, $\lambda_{NP\rightarrow P} = " +str(p_arrive_i) + "$"
                    ,ylabel=ylab
                    )

    # x axis label
    the_fig.fig.text(
            x=0.52
            ,y = 0.02
            ,s="Time, $t$"
            ,fontsize=18
            ,horizontalalignment="center"
            ,transform = the_fig.fig.transFigure)
    
    # y axis label
    the_fig.fig.text(
            x=0.02
            ,y = 0.52
            ,s="Hormone level, $z$"
            ,rotation=90
            ,fontsize=18
            ,horizontalalignment="center"
            ,transform = the_fig.fig.transFigure)

    the_fig.close(tight=True)

zmax = 1
dmax = 1
aP = [ 0.5, 2]
ad = [ 0.5, 2]


# damage decay
r = [1]
#u = [0.25, 0.5, 0.75,1]
u = [ 0.9 ]

filename = "stress_iteration_overview"

for aP_i in aP:
    for ad_i in ad:
        for r_i in r:
            for u_i in u:

                filename_sub = filename +\
                        "_ad_" + str(ad_i) +\
                        "_aP_" + str(aP_i) +\
                        "_r_" + str(r_i) +\
                        "_u_" + str(u_i)

                stress_multipanel(
                        {"init_stress_influx": 0
                            ,"init_feedback": 0
                            ,"init_influx": 0
                            ,"cue_P": 0
                            ,"dmax": dmax
                            ,"zmax":zmax
                            ,"ad":ad_i
                            ,"aP":aP_i
                            ,"u":u_i
                            ,"r":r_i
                            }
                        ,filename_sub + "_inits0.pdf")
