#!/usr/bin/env Rscript
# get rid of everything there
rm(list=ls())

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("filesstrings"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("viridis"))

# some general functions
source("stress_file_functions.r")


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
        pattern="stress"
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

# first plot: hormone
p1 <- ggplot(data.attack.acute) +
        geom_line(aes(y=hormone, x=time)) +
        ggtitle("Acute stress")

# second plot: hormone
p2 <- ggplot(data.attack.acute) +
        geom_line(aes(y=damage, x=time))

# first plot: hormone
p3 <- ggplot(data.attack.chronic) +
        geom_line(aes(y=hormone, x=time)) +
        ggtitle("Chronic stress")

# second plot: hormone
p4 <- ggplot(data.attack.chronic) +
        geom_line(aes(y=damage, x=time))


        #test.u <- as_tibble(expand.grid(t = seq(1,50,1)
        #        ,d = seq(0,1000,1)))
        #
        #test.u$hormone <- rnorm(n=nrow(test.u))
        #print(nrow(test.u))


# now a levelplot showing the stress levels
p5 <- ggplot(data=data.strategy) +
        geom_tile(aes(x = t, y = d, fill = hormone)) +
        scale_x_continuous(limits = c(0,50), expand=c(0,0)) +
        scale_y_continuous(limits = c(0,1000), expand=c(0,0)) +
        scale_fill_viridis() +
        ggtitle("Hormone strategy over tau and damage") 


# make annotation
param.str <- paste0(
        "pLeave: ",params["pLeave"],",",
        " pArrive: ",params["pArrive"],",",
        " repair: ",params["repair"],",",
        " K: ",params["K"])

(p1 / p2 / p3 / p4 / p5) +
    plot_annotation(param.str)

ggsave(filename=plot.file.name)
