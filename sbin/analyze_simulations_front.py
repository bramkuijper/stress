#!/usr/bin/env python

# to analyze individual-based sims, but then parameters are given upfront

import os, re, sys

first = True


def analyze_parameters(lines,first=False):

    pars = {}

    for line in lines:

        splitpars = line.split(";")

        if len(splitpars) > 1:
            pars[splitpars[0]] = splitpars[1]

    return(pars)

def analyze_data(lines):

    data = [ [ float(celli) ] for celli in lines[0].split(";")[0:-1] ]

    # loop through the lines to collect the data
    for line in lines[1:]:
        splitline = line.split(";")[0:-1]

        for i in range(0,len(splitline)):
            data[i].append(float(splitline[i]))

    # now take averages

    avgs = []
    for i in range(0,len(data)):
        avgs.append(sum(data[i])/len(data[i]))

    return(avgs)

def analyze_file(filename):

    global first;

    # open file; read data
    f = open(filename)
    fl = f.readlines()
    f.close

    parameter_endline = 0

    # first see until what line the parameters stretch
    for idx, line_i in enumerate(fl):
        if re.match("^generation",line_i) != None:
            parameter_endline = idx
            break

    # something went wrong. Could not find 
    # the header. return
    if parameter_endline == len(fl):
        return

    data_endline = -1

    linelist = list(range(1,len(fl)))
    linelist.reverse()

    # now find the last dataline at the end of the file
    for line_number in linelist:

        # found a line with data, this is my last line
        if re.match("^\d.*",fl[line_number]) is not None:
            data_endline = line_number;
            break

    flhead = fl[parameter_endline]

    parameters = analyze_parameters(fl[0:parameter_endline])

    if first:
        print(";".join(parameters.keys()) + ";" + flhead.strip() + "file")
        first = False

    print(";".join(parameters.values()) + ";" + fl[data_endline].strip() + filename)


for root, dirs, files in os.walk(sys.argv[1]):

    for name in files:
        if re.match("(sim|iter).*\d$",name) != None:
            data = analyze_file(os.path.join(root,name))


