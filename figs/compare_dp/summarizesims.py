#!/usr/bin/env python3

import os
import re
import sys
import pandas as pd

# class to summarize simulation output
# where csv files with data + parameters are parsed
# so that a summary file is generated with the last
# line of the data + all parameters + filename for each file
class SummarizeSims:
    def __init__(self
                 ,path="."
                 ,pattern=r"(sim.*\d+$|iter.*\d+)"
                 ,recursive=False
                 ,sep=";"
                 ,parameters_first=False):
        """
        Initializes the summarize sims class

        Parameters
        ----------
        path : str, optional
            The path in which files need to found. The default is ".".
        pattern : str, optional
            Regular expression matching the filenames that need to be 
            processed. The default is r"sim|iter".
        recursive : boolean, optional
            Whether the directory should be searched recursively. The
            default is False
        sep : str, optional
            The separator used in the file. The default is ";"
        parameters_first : boolean, optional
            Whether the listing of the parameters occurs 
            at the start of the file

        Returns
        -------
        None.

        """
        
        self.first = True
        self.init = False
        self.path = path
        self.num_semicol_header = 0
        self.pattern = pattern
        self.recursive = recursive
        self.sep = sep
        

        self.full_data = None
        
        for root, dirname, files in os.walk(self.path):
            
            for file in files:
                if re.search(self.pattern,file) is not None:
#                    print(os.path.join(root,file))
                    self.analyze_file(filename=os.path.join(root,file))

        self.output_full_data()

    def output_full_data(self):

        # write the data to stdout
        if self.full_data.shape[0] > 0:
            print(self.full_data.to_csv(
                    path_or_buf=None
                    ,sep=self.sep
                    ,index=False))

    # analyze parameters at the end of the file
    def analyze_parameters_old(self, lines):
        """
        Runs through the parameter part and analyzes it

        Parameters
        ----------
        lines : list
            list of strings, where each string is a line from the original
            simulation file.
        
        Returns
        -------
        None.

        """
    
        pars = {}
    
        for line in lines:
            mobj = line.split(self.sep)
    
            if len(mobj) > 1:
                pars[mobj[0]] = mobj[1]
    
        return(pars)

    def analyze_data(
            self
            ,filename
            ,linenumber_start=0
            ,linenumber_end=None):

        the_data = pd.read_csv(
                filepath_or_buffer=filename
                ,skiprows=linenumber_start
                ,nrows=linenumber_end
                ,sep=self.sep)

        return(the_data.iloc[[-1]])
    
    def analyze_parameters_new(
            self 
            ,linenumber_start
            ,filename
            ,linenumber_end=None):

        assert(type(linenumber_start) == type(30))
        assert(type(filename) == type("adf"))
        assert(linenumber_end > linenumber_start)

        # read in the parameter data as if it were
        # a csv file
        param_data = pd.read_csv(
                filepath_or_buffer=filename
                ,sep=self.sep
                ,skiprows=linenumber_start
                ,nrows=linenumber_end - 1
                ,header=None
                ,dtype=str
                ,usecols=[0,1]).T

        param_data.columns = param_data.iloc[0]
        param_data.drop(index=0,axis=0, inplace=True)

        return(param_data) 
    
    # processes the first line headers
    # when making line headers for initial values
    def process_first_line(self, line):
    
        # get the column names and split them into a list
        line_cols = line.strip().split(self.sep)
    
        new_cols = ""
    
        for colname in line_cols:
            if not colname:
                continue
    
            new_cols += colname + "_t_0" + self.sep
    
        return(new_cols)

    def find_last_data_line(self, filename):

        with open(filename) as infile:
            linelist = infile.readlines()

            # make a reverse range and loop through lines
            linerange = range(len(linelist),0,-1)

            for line_idx in linerange:
                if re.search(r"^\d",linelist[line_idx]) is not None:
                    return(line_idx)


    
    def analyze_file(self,filename, sep=";"):

        # two options (ignoring empty lines)
        # 1. header line, data, then parameters
        # 2. parameters, header line, data 

        # parameters recognizable as they contain
        # non-numerical data at start. Hence multiple
        # non-numerical starts consecutively means parameters

        non_data_lines_start = []
        data_start = None

        # open the file
        with open(filename) as infile:

            # loop through lines from start
            for line_idx, line in enumerate(infile):
                # skip empty lines
                if line in ["\n","\r\n"]:
                    continue

                # check whether line contains parameters
                # or column headings
                if re.search(r"^\s*\d",line) is None:
                    non_data_lines_start.append(line_idx)
                else:
                    # ok we hit a line containing numerical
                    # data
                    # we can break as we know enough
                    data_start = line_idx
                    break
            else:
                # ok no data was found, as a break should
                # have occurred somewhere
                # we will ignore this file
                return None


        # we are done sniffing the data, we can now make a choice
        # if more than one non_data_line, we got a parameter listing
        # at the start
        if len(non_data_lines_start) > 1:

            # parameter listing at the start, call analyze parameters function
            parameters = self.analyze_parameters_new(
                    linenumber_start=non_data_lines_start[0]
                    ,linenumber_end=non_data_lines_start[-2] # not the last element as that is the data header
                    ,filename = filename)

            data = self.analyze_data(
                    filename=filename
                    ,linenumber_start=non_data_lines_start[-1]
                    )

            # combine both parameters and data
            data_params_combined = pd.concat(
                    [parameters.reset_index(drop=True)
                        ,data.reset_index(drop=True)]
                    ,axis=1)
            
            # add current filename
            data_params_combined["file"] = filename
            

            if self.full_data is None:
                self.full_data = data_params_combined
            else:
                self.full_data = self.full_data.append(
                        other=data_params_combined
                        ,ignore_index=True
                        ,sort=False
                        )

        else: # parameters at the end apparently

            # get a list of all non-data-lines at the end of the file
            last_data_line = self.last_data_line(filename)

            # collect the data
            data = self.analyze_data(
                    linenumber_start=non_data_lines_start[-1]
                    ,linenumber_end=last_data_line - 1
                    ,filename=filename)
               
            parameters = self.analyze_parameters_new(
                    linenumber_start=last_data_line + 1
                    ,linenumber_end=None
                    ,filename = filename)

#        # indicator variable whether we are at first line
#        # of the file to be read
#        firstline = True
#    
#        flhead = ""
#    
#        lc = 0
#    
#        # indicator variable whether we are at the part
#        # involving parameters
#        parameter_part = False
#    
#        parameter_lines = []
#    
#        # the line where the parameter output
#        # starts
#        parline = -1
#    
#        # store the last line of data
#        last_data_line = ""
#    
#        # store the first line of data
#        # in case we need initial values too
#        first_data_line =""
#    
#        # the header of the resulting data file
#        flhead = ""
#    
#        # open the file and read the stuff
#        with open(filename) as infile:
#            for line in infile:
#    
#                # see whether we also have to store the initial values
#                if self.init:
#                    if firstline:
#                        flhead += process_first_line(line, sep=sep)
#    
#                # update line count
#                lc += 1
#    
#                # get the first line of data
#                if lc == 2:
#                    first_data_line = line.strip()
#    
#                # if this is the first line store
#                # the header
#                if firstline:
#                    flhead += line.strip()
#                    firstline = False
#    
#                # if this is any other line starting
#                # with a numerical value store the line
#                # as it might be potentially the last one
#                elif re.match("^\d",line):
#                    last_data_line = line
#    
#                # hold this until we have the parameter file
#                if not parameter_part:
#                    if re.match("^\n",line) is not None:
#                        parline = lc
#                        parameter_part = True
#                        parameter_lines += [line.strip()]
#                elif parameter_part:
#                    parameter_lines += [line.strip()]
#    
#        if parline < 1:
#            return
#    
#        parameters = self.analyze_parameters_new(
#                linenumber_start=parline
#                ,filename=filename)
#    
#        # prepare the data to be printed
#        data = ""
#    
#        if self.init:
#    
#            # do error checking in terms of the csv values
#            count_semicol = len(re.findall(sep, first_data_line))
#            count_semicol_data = len(re.findall(sep, last_data_line))
#    
#            assert(count_semicol == count_semicol_data)
#    
#            data += first_data_line.strip()
#    
#        data += last_data_line.strip()
#    
#    
#        if first:
#    
#            flhead = flhead.strip()
#    
#            # check if flhead has trailing comma
#            if re.search(sep.encode("unicode-escape").decode() + "$"
#                    ,flhead) is None:
#                flhead = flhead + sep
#    
#            header_line = sep.join(parameters.keys()) + sep + flhead + "file"
#    
#            # count number of occurrences of semicolon for error checking
#            num_semicol_header = len(re.findall(sep,header_line))
#    
#            print(header_line)
#    
#            first = False
#    
#        # check if data has trailing separator
#        # if not we'll have to add one
#        if re.search(sep.encode("unicode-escape").decode() + "$"
#                ,data) is None:
#            data = data + sep
#    
#        data_line = sep.join(parameters.values()) + sep + data +  filename
#    
#        # count number of occurrences of semicolon for error checking
#        num_semicol_data = len(re.findall(sep,data_line))
#    
#        assert(num_semicol_header == num_semicol_data)
#        print(data_line)
