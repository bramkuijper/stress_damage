#!/usr/bin/env bash

find . -iname "stress_output*.txt" -print0 | xargs -0 -P8 -I% ./single_plot.r %
