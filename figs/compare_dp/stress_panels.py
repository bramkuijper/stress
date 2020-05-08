#!/usr/bin/env python3

# plots the first generation's outcome
from pythontools import multipanel
import pandas as pd
from matplotlib import cm
import os.path
import sys
import re

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

def get_dp_iter(filename):
    
    line_number = None
    

    
    with open(filename) as dpfile:
        
        lines = dpfile.readlines()
           
        linerange = range(len(lines)-1,0,-1)
        
        for line_idx in linerange:
            
            if re.search(pattern="^\d", string=lines[line_idx].strip()) is not None:
                line_number = line_idx - 1
                break
    
        else:
            return None
    
 
    # get the contents of the file
    dp_iter_dat = pd.read_csv(
        filepath_or_buffer=filename
        ,nrows=line_number
        ,sep="\t")
    
      
    return dp_iter_dat

# from a summary data frame obtain filenames
# and extract iteration data from those
def get_iter_data(df, iter_dir="."):
    
    assert(os.path.exists(iter_dir))
    
    assert(df.shape[0]>0)

    # get list of filenames
    file_list = list(df["file"].unique())
    
    assert(len(file_list) > 0)

    iter_df = None


    # now open files and extract data
    # for each combination of parameters there are likely
    # multiple 
    for replicate, file_name in enumerate(file_list):

        # get the filename of the filename containing
        # the iterations
        file_base_name = os.path.basename(file_name)

        iter_file_name = os.path.join(iter_dir, file_base_name + "iters.csv")

        try:
            # read the iteration data
            subdf = pd.read_csv(iter_file_name,sep=";")

            subdf["replicate"] = replicate

            if iter_df is None:
                iter_df = subdf
            else:
                iter_df = iter_df.append(subdf, ignore_index=True)

        except:
            print("meh file " + iter_file_name + " won't exist")
            
    return(iter_df)



#### function declarations ####

def single_panel(
        data_dp
        ,data_sims
        ,multipanel_obj
        ,row_i=0
        ,col_i=0
        ,min_time_iter = 0
        ,max_time_iter = 200
        ,newzero = 0
        ,ylim=[-0.05,1.05]
        ,xlim=[-30,200]
        ,xticks=True
        ,xlab=""
        ,ylab=""
        ,title=""
        ):
    
    
    """
    Plot a single plotting panel

    Parameters
    ----------
    data_dp : pandas dataframe
        The applicable DP summary data contained in a dataframe.
    data_sims : pandas dataframe
        DESCRIPTION.
    row_i : int, optional
        DESCRIPTION. The default is 0.
    col_i : int, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    None.

    """
    sim_line_color = "#46a9ff"
    sim_line_alpha = 0.2
    dp_line_color = "black"
    dp_line_alpha = 0.5
#    thin_line_color = "grey"
    
    
    ######## make a single panel ########
    the_axis = multipanel_obj.start_block(
                row=row_i
                ,col=col_i)
    
    
    # get colors
    the_color_map = cm.get_cmap("tab10")
    
    iteration_data = None
        
    if data_sims.shape[0] > 0:
        dirname = os.path.dirname(str(data_sims["file"].values[0]))
    
        # obtain the iteration file data
        iteration_data = get_iter_data(
            data_sims
            ,iter_dir=dirname)

    if type(iteration_data) == type(pd.DataFrame()):

        iteration_data = iteration_data.loc[
                (iteration_data["time"] <= max_time_iter)
                &
                (iteration_data["time"] >= min_time_iter)
                    ].copy(deep=True)

        iteration_data["time"] = iteration_data["time"] + newzero

        # get list of unique replicates per parameter combination
        replicates = list(iteration_data["replicate"].unique())
       
        # get list of unique individuals per replicate
        individuals = list(iteration_data["individual"].unique())

        print("number replicates per plot: " + str(len(replicates)))

        for replicate_i in replicates:
            for individual_i in individuals:

                # select the subset for individual i and replicate j
                iteration_subset = iteration_data.loc[
                        (iteration_data["replicate"] == replicate_i)
                        &  (iteration_data["individual"] == individual_i)].copy(deep=True)
                
                # sort the subset
                iteration_subset.sort_values(by="time",inplace=True)

                the_axis.plot(iteration_subset["time"]
                        ,iteration_subset["hormone"]
                        ,color=sim_line_color
                        ,linewidth=0.5
                        ,alpha=sim_line_alpha
                        ,label="_nolabel")

    print("dp data SHAPE: " + str(data_dp.shape[0]))
                
    if data_dp.shape[0] > 1:
        data_dp = data_dp.iloc[[1],:]

    print("dp data SHAPE now: " + str(data_dp.shape[0]))
                
    if data_dp.shape[0] == 1:
                    
        dp_iter = get_dp_iter(
            filename=str(data_dp["file"].values[0]))

        if "t" not in dp_iter.columns.values:
            print(data_dp["file"].values[0])
            print(dp_iter.columns.values)
            raise
        
        dp_iter = dp_iter.loc[(dp_iter["t"] >= xlim[0])
                              & (dp_iter["t"] <= xlim[1])]
        
     
        timepoints = [min_time_iter + newzero, -1, 0] + list(dp_iter["t"].values[:])
        hormone = [dp_iter["hormone"].values[-1]
                   ,dp_iter["hormone"].values[-1]
                   ,dp_iter["hormone"].values[0]
                   ] + list(dp_iter["hormone"].values[:])
        
        hormone = [ i/1000 for i in hormone ]
        
        # now plot the dp data
        the_axis.plot(timepoints
                      ,hormone
                      ,color=dp_line_color
                      ,alpha=dp_line_alpha
                      ,linewidth=4)          

    ylab = ""

    # end the figure
    multipanel_obj.end_block(
            the_axis
            ,ylim=ylim
            ,xlim=xlim
            ,y_ticks_minor = 1
            ,x_ticks_minor = 1
            ,x_ticks_major_multiple = 25
            ,y_ticks_major_multiple = 0.25
            ,xticks=xticks
            ,yticks=col_i == 0
            ,title=title
            ,xlabel=xlab
            ,ylabel=ylab
            ,loc_title=True
            ,loc_title_pos=[-0.025,1.05]
            )




# simple function to use a dict to query a dataframe
def make_query(the_dict):

    qry_list = []

    for key, value in the_dict.items():
        qry_list += [str(key) + " == '" + str(value) + "'"]

    qry = " & ".join(qry_list)

    return(qry)

def check_columns(
        param_dict
        ,data
        ):
    
    colnames = set(data.columns.values)
    param_colnames = set(param_dict.keys())
    
    if len(param_colnames - colnames) > 0:
        raise Exception("Some column names in query do not exist in the dataset " + str(param_colnames - colnames))


######## general figure function ########




def stress_multipanel(
        param_array
        ,sim_data
        ,dp_data
        ,filename
        ,title_array=None
        ,min_time_iter = 0
        ,max_time_iter = 150
        ,newzero = 0
        ,xlim=[0,100]
        ):
    
    nrows, ncols = param_array.shape
    
    rows = [1 for i in range(nrows)]
    cols = [1 for i in range(ncols)]
    
    
    # use the dictionary to make a subset of the 
   
    # initialize multipanel fig
    the_fig = multipanel.MultiPanel(
            panel_widths=cols
            ,panel_heights=rows
            ,filename=filename
            ,width=10
            ,wspace=0.2
            ,hspace=0.3
            )
    
    try:
        
        # now loop over rows and columns
        for row_i, row in enumerate(param_array):
            for col_i, single_panel_dictionary in enumerate(row):
                
                # check whether current item is indeed a dictionary
                assert(type(single_panel_dictionary) == type({}))
                
                check_columns(
                    param_dict=single_panel_dictionary
                    ,data=sim_data)
                
                check_columns(
                    param_dict=single_panel_dictionary
                    ,data=dp_data)
                        
                # make a query string
                query_string = make_query(single_panel_dictionary)
            
                sim_subset = sim_data.query(query_string)
                dp_subset = dp_data.query(query_string)
                
                title = ""
                
                if title is not None:
                    
                    title= title_array[row_i][col_i]
                    

                single_panel(
                    data_dp=dp_subset
                    ,data_sims=sim_subset
                    ,multipanel_obj=the_fig
                    ,row_i = row_i
                    ,col_i = col_i
                    ,min_time_iter=min_time_iter
                    ,max_time_iter=max_time_iter
                    ,newzero=newzero
                    ,xlim=xlim
                    ,title=title
                    )

                # check whether there is iteration data for this
                # data set
    
                # if yes, plot it
    except:
        print("something went wrong...")
        the_fig.close()
        raise

    # x axis label
    the_fig.fig.text(
            x=0.52
            ,y = 0.02
            ,s="Time, $t$"
            ,fontsize=18
            ,horizontalalignment="center"
            ,transform = the_fig.fig.transFigure)

    # make string of parameter values and put on plot
    translate_parameters = {
            "init_feedback":r"$\mathrm{Feedback}_{t_{0}}$:"
            ,"ad":r"$a_{d}$:"
            ,"aP":r"$a_{P}$:"
            ,"u":r"$u$:"
            ,"r":r"$r$:"
            }

    # y axis label
    the_fig.fig.text(
            x=0.05
            ,y = 0.52
            ,s="Hormone level, $z$"
            ,rotation=90
            ,fontsize=18
            ,horizontalalignment="center"
            ,verticalalignment="center"
            ,transform = the_fig.fig.transFigure)

    print("close")
    the_fig.close(tight=True)
