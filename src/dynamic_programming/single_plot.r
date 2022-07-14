#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("filesstrings"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library("here"))

# load file with general functions
script.dir <- here()
source(file.path(script.dir,
                "src/dynamic_programming/stress_file_functions.r"))


#### WORK OUT THE FILENAMES ####

# get the file containing the data which we want
# to plot from the command line
args = commandArgs(trailingOnly=TRUE)
strategy.file.name <- args[[1]]

# obtain file name without final extension (.txt)
# using files::before_last_dot()
strategy.file.base <- before_last_dot(strategy.file.name) 

# from given file name obtain the simAttack file
# by replacing stress for simAttacks
attack.file.name <- gsub(
        pattern="stress_strategy"
        ,replacement="simAttacks"
        ,x=strategy.file.name)

plot.file.name <- paste0("plot_",basename(strategy.file.base),".pdf")



### GET THE ACTUAL DATA ###

# get the parameters of this run
params <- read.parameters(strategy.file.name)

# get the strategy data
data.strategy <- read.stress.file(strategy.file.name)

# get the simAttack data
data.attack <- read.attack.file(attack.file.name)

data.attack.acute <- filter(data.attack, type == "acute")
data.attack.chronic <- filter(data.attack, type == "chronic")


### PLOTTING ###
list.plots <- list()

# first plot: hormone
list.plots <- c(list.plots,list(ggplot(data.attack.acute) +
        geom_line(aes(y=hormone, x=time)) +
        theme_classic() +
        ggtitle("Acute stress")))

# second plot: damage
list.plots <- c(list.plots, list(ggplot(data.attack.acute) +
        theme_classic() +
        geom_line(aes(y=damage, x=time))))

# first plot: hormone
if (nrow(data.attack.chronic) > 0)
{
    list.plots <- c(list.plots, list(ggplot(data.attack.chronic) +
        geom_line(aes(y=hormone, x=time)) +
        theme_classic() +
        ggtitle("Chronic stress")))

    # second plot: hormone
    list.plots <- c(list.plots, list(ggplot(data.attack.chronic) +
            theme_classic() +
            geom_line(aes(y=damage, x=time))))
}


# now a levelplot showing the stress levels
list.plots <- c(list.plots, list(ggplot(data=data.strategy) +
        geom_tile(aes(x = t, y = d, fill = hormone)) +
        scale_x_continuous(limits = c(-0.5,50.5), expand=c(0,0)) +
        scale_y_continuous(limits = c(0,1000), expand=c(0,0)) +
        scale_fill_viridis(option="magma") +
        xlab("Time since last attack, tau") +
        theme_classic() +
        ggtitle("Hormone strategy over tau and damage")))

# make annotation
param.str <- paste0(
        "pLeave: ",params["pLeave"],",",
        " pArrive: ",params["pArrive"],",",
        " repair: ",params["repair"],",",
        " pAtt: ",params["pAttack"],",",
        " maxTs: ",params["maxTs"],","
        )

autocorr <- 1.0 - (params["pLeave"] + params["pArrive"])
risk <- params["pArrive"] / (params["pLeave"] + params["pArrive"])

param.str <- paste0(param.str,
        " risk: ",risk,",",
        " autocorr: ",autocorr,",")

if ("n_repro_bout" %in% names(params))
{
    param.str <- paste0(param.str,
            " bout: ",params["n_repro_bout"])
}



# find fecundity values
Kcolumns <- grepl("^K",names(params))
fec.cols <- names(params)[Kcolumns]

for (col in fec.cols) 
{
    param.str <- paste0(param.str,", ",col,": ",params[col])
}


all.plots <- wrap_plots(list.plots, ncol=1) + 
    plot_annotation(title=param.str,
            caption=strategy.file.name,
            theme=theme(text=element_text(size=8)
                    )) 
#    plot_layout(heights = unit(c(rep(1,times=4),8),rep("null",times=5)))

ggsave(filename=plot.file.name,height=12)
