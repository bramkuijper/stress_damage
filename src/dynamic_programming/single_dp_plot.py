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
import numpy as np
mpl.use("pgf")


#### fonts etc

#fontpath = "/System/Library/Fonts/Supplemental/" 
#fontpath = "/System/Library/Fonts/Supplemental/" 
#
#pgf_with_custom_preamble = {
#    "font.family": "serif", # use serif/main font for text elements
#    "text.usetex": True,    # use inline math for ticks
#    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
#    "pgf.preamble": "\n".join([
#        r"\usepackage{units}",         # load additional packages
#        r"\usepackage{mathspec}",         # load additional packages
#        r"\setmainfont[" +\
#            "Path = " + fontpath + "," +\
#            "UprightFont = * ," +\
#            "ItalicFont = *Italic ," +\
#            "BoldFont = *Bol," +\
#            "Extension = .otf]{STIXGeneral}",
#        r"\setmathsfont(Digits,Latin,Greek)[" +\
#            "Path = " + fontpath + "," +\
#            "UprightFont = * ," +\
#            "ItalicFont = *Italic ," +\
#            "BoldFont = *Bol," +\
#            "Extension = .otf]{STIXGeneral}",
#        r"\setmathrm[" +\
#            "Path = " + fontpath + "," +\
#            "UprightFont = * ," +\
#            "ItalicFont = *Italic ," +\
#            "BoldFont = *Bol," +\
#            "Extension = .otf]{STIXGeneral}",
#         ])
#}

#mpl.rcParams.update(pgf_with_custom_preamble)

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

############ pivot tables of the histograms ############

def generate_pivot(the_data,x, y, z):
    the_pivot = the_data.pivot_table(values=z, index=y, columns=x)

    x, y = np.meshgrid(the_pivot.columns.values, the_pivot.index.values)

    z = the_pivot.values

    return(x, y, z)



#### get the stress_data file name
file_name = sys.argv[1]

# do some checking
if re.search(r"stressL.*txt$",file_name) == None:
    raise Exception("Incorrect filename, give me one of the stressLXXX.txt ones.")

# get the footer of the dataset
end_line, pardict = find_footer(file_name)

maxH = pardict["maxH"]
maxD = pardict["maxD"]

seasonal = False

if "maxTs" in pardict.keys():
    maxTs = pardict["maxTs"]

    seasonal = True

    tsval = [ 0, round(maxTs/4), round(maxTs/2), round(maxTs * 0.75), round(maxTs) ]




##### read in the data #####
skiprows=2

stress_data = pd.read_csv(filepath_or_buffer=file_name
        ,sep="\t"
        ,skiprows=skiprows
        ,nrows=end_line-skiprows)

# get all the damage values
damage_vals = list(stress_data["d"].unique())
damage_vals.sort()

damage_filter = [0,0.25, 0.5,0.75,1]

damage_val_select = []
for damage_filter_i in damage_filter:
    damage_val_select += [damage_vals[round(damage_filter_i*(len(damage_vals) - 1))]]
print(damage_val_select)

# read in the simulated attack data
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

# at least 2 rows
# for the simulation
nrows = 2

# if this is a simulation with seasonality
# we will plot as many imshow plots()
# as there are selected levels of damage 
# otherwise we just produce a single line plot
nrows += len(damage_val_select) if seasonal else 1


the_fig = multipanel.MultiPanel(
        panel_widths=[1]
        ,panel_heights=[1 for x in range(0,nrows)]
        ,filename="plot_" + file_root + ".pdf"
        ,hspace=0.3
        ,width=8
        ,height=15
        )

rowctr = 0

if seasonal:

    for d_i in damage_val_select:

        the_axis = the_fig.start_block(
                row=rowctr
                ,col=0)

        stress_subset = stress_data[stress_data["d"] == d_i].copy(deep=True)

        # make pivot table to fit 
        # the imshow or contourplot
        (x, y, z) = generate_pivot(stress_subset
                ,x="t"
                ,y="ts"
                ,z="hormone")

        the_axis.imshow(z,
                    extent=[x.min(),
                                x.max(),
                                y.min(),
                                y.max()],
                    origin="lower",
                    aspect="auto",
                    interpolation="none",
                    cmap=cm.get_cmap(name="viridis"))

        the_fig.end_block(ax=the_axis
                ,ylabel="Ts"
                ,xticks=damage_val_select.index(d_i) == len(damage_vals)
                ,yticks=True
                ,xlabel="Time"
                ,title=r"Damage = " + str(d_i))

        rowctr += 1

else:
    the_axis = the_fig.start_block(
            row=rowctr
            ,col=0)

    for d_i in damage_val_select:
        
        d_sub = stress_data[stress_data["d"] == d_i]

        the_axis.plot(d_sub["t"]
                ,d_sub["hormone"]
                ,label=r"$d = " + str(d_i) + "$"
                )

    the_axis.legend()

    the_fig.end_block(ax=the_axis
            ,ylabel="Hormone"
            ,xticks=True
            ,yticks=True
            ,xlabel="Time")

    rowctr += 1



## plot simulated attacks:hormone
the_axis = the_fig.start_block(
        row=rowctr
        ,col=0)

rowctr += 1

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
        ,yticks=True
        ,ylim=[-0.01*maxH,maxH+0.01*maxH])

## plot simulated attacks: damage
the_axis = the_fig.start_block(
        row=rowctr
        ,col=0)

rowctr += 1

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
        ,yticks=True
        ,ylim=[-1,maxD+1])


the_fig.close(tight=True)
