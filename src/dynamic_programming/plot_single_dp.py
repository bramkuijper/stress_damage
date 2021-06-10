#!/usr/bin/env python3

from pythontools import multipanel
import itertools
import pandas as pd
import sys
import re
import os.path
import argparse
from matplotlib import cm
import matplotlib.lines as lines
import matplotlib.patches as mpatches
import matplotlib as mpl
mpl.use("pgf")


#### fonts etc

fontpath = "/System/Library/Fonts/Supplemental/" 

pgf_with_custom_preamble = {
    "font.family": "serif", # use serif/main font for text elements
    "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": "\n".join([
        r"\usepackage{units}",         # load additional packages
        r"\usepackage{mathspec}",         # load additional packages
        r"\setmainfont[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
        r"\setmathsfont(Digits,Latin,Greek)[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
        r"\setmathrm[" +\
            "Path = " + fontpath + "," +\
            "UprightFont = * ," +\
            "ItalicFont = *Italic ," +\
            "BoldFont = *Bol," +\
            "Extension = .otf]{STIXGeneral}",
         ])
}

mpl.rcParams.update(pgf_with_custom_preamble)

def process_params(line_array):

    pardict = {}

    for line in line_array:
        lst = re.split(pattern=":",string=line)

        if len(lst) > 1:
            pardict[lst[0]] = float(lst[1].strip())

    return(pardict)

# find the line where the stress_data ends
# and return a dict with parameters
def find_footer(file_name):

    lines = None

    with open(file_name) as lfile:
        lines = lfile.readlines();

    # loop over the *reversed* line range
    linerange = range(len(lines)-1,0,-1)

    found_params = False

    param_dict = {}

    for line_idx in linerange:

        if not found_params and re.search(pattern="PARAMETER", string=lines[line_idx]) != None:
            param_dict = process_params(lines[line_idx+1:])
            found_params = True

        if found_params and re.search(r"^nIterations",lines[line_idx]) != None:
            return (line_idx - 2, param_dict)

# end find_footer(file_name):

# remove trailing tabs from some files
# as pd.read_csv chokes on it
def rm_trailing_tabs(file_name):

    with open(file=file_name) as the_file:
        fl = the_file.read()

    fl = re.sub(pattern=r"\t$"
            ,repl=""
            ,string=fl
            ,flags=re.MULTILINE)

    with open(file=file_name, mode="w") as the_file:
        the_file.write(fl)


#### get the stress_data
file_name = sys.argv[1]

# do some checking
if re.search(r"stressL.*txt$",file_name) == None:
    raise Exception("Incorrect filename, give me one of the stressLXXX.txt ones.")

# get the footer of the dataset
end_line, pardict = find_footer(file_name)

maxH = pardict["maxH"]
maxD = pardict["maxD"]

skiprows=2

stress_data = pd.read_csv(filepath_or_buffer=file_name
        ,sep="\t"
        ,skiprows=skiprows
        ,nrows=end_line-skiprows)

file_name_sim_attack = re.sub(
        pattern="stressL"
        ,repl="simAttacksL"
        ,string=file_name)

rm_trailing_tabs(file_name_sim_attack)

sim_attack_data = pd.read_csv(
        filepath_or_buffer=file_name_sim_attack
        ,sep="\t")

# start the figure
# file name without extension
path, file_name_last = os.path.split(file_name)

file_root, file_ext = os.path.splitext(file_name_last)

the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1,1,1,1]
        ,filename="plot_" + file_root + ".pdf"
        ,hspace=0.3
        ,width=8
        ,height=5
        )
    
the_axis = the_fig.start_block(
        row=0
        ,col=0)

the_axis.plot(stress_data["t"]
        ,stress_data["d"])

the_fig.end_block(ax=the_axis
        ,ylabel="Damage"
        ,ylim=[0,maxD])

## plot hormone over time
the_axis = the_fig.start_block(
        row=1
        ,col=0)


the_axis.plot(stress_data["t"]
        ,stress_data["hormone"])

the_fig.end_block(ax=the_axis
        ,ylabel="Hormone"
        ,ylim=[0,maxH]
        ,xticks=True
        ,xlabel="Time")




## plot simulated attacks:hormone
the_axis = the_fig.start_block(
        row=2
        ,col=0)

# plot hormone levels
the_axis.plot(sim_attack_data["time"]
        ,sim_attack_data["hormone"]
        ,color="blue")

# now plot the timepoints at which an attack took place
attack_data_fltr = sim_attack_data[sim_attack_data["attack"] == 1]

nrow_attack = attack_data_fltr.shape[0]

if nrow_attack > 0:
    the_axis.plot(attack_data_fltr["time"],
            [1.0/10 * maxH for x in range(0, nrow_attack,1)],
            linestyle="",
            marker="o",
            markersize=3,
            markerfacecolor="red",
            markeredgecolor="red"
            )

the_fig.end_block(ax=the_axis
        ,ylabel="Hormone"
        ,ylim=[0,maxH])

## plot simulated attacks: damage
the_axis = the_fig.start_block(
        row=3
        ,col=0)

the_axis.plot(sim_attack_data["time"]
        ,sim_attack_data["damage"]
        ,color="magenta")

if nrow_attack > 0:
    the_axis.plot(attack_data_fltr["time"],
            [2 for x in range(0, nrow_attack,1)],
            linestyle="",
            marker="o",
            markersize=3,
            markerfacecolor="red",
            markeredgecolor="red"
            )

the_fig.end_block(ax=the_axis
        ,ylabel="Damage"
        ,xticks=True
        ,xlabel="Time"
        ,ylim=[0,maxD])


the_fig.close(tight=True)
