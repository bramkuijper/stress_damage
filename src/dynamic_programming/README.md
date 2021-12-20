# Evolution of stress responses in the face of damage and when reproduction is seasonal

In Taborsky et al (2021) Trends Ecol Evol, we investigated the evolution of stress responses in an autocorrelated environment without the accumulation of stress-related damage.

Here we relax this assumption and allow damage `d` to accumulate, where damage accumulation increases with hormone expression `h` and decreases with repair `r`.

## Compiling the code
The source code of the dynamic programming model is written in `C++` and is available in the file `stress_damage_lh.cpp`. Compiling the programme can be done by running
```
make
```
which yields the executable `stress_damage_lh.exe`.

## Running the code
The programme expects command-line parameters (for more info on the command line params, see the function `init_params()` in the c++ file `stress_damage_lh.cpp`. These parameters could be provided in a plain-text bash shell file, such as the example file provided in `example_run.sh`. In a bash shell, one can run:
```
bash run_single_dp.sh
```

## Algorithm to vary the parameters
Parameters can also be varied by using a bunch of `for`-loops provided in the Python3 script `vary_params.py`. A file called `runs.sh` containing a series of commands to run replicate simulations can be obtained by running
```
./vary_params.py > runs.sh
```

The file will look like
```
./stress_damage_lh.exe 0.95 0.05 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.9 0.1 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.8 0.2 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.855 0.045 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.81 0.09 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.72 0.18 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.665 0.035 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.63 0.07 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.56 0.14 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.475 0.025 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.45 0.05 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.4 0.1 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.285 0.015 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.27 0.03 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.24 0.06 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.095 0.005 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.09 0.01 0.5 1.0 0.001 0.01
./stress_damage_lh.exe 0.08 0.02 0.5 1.0 0.001 0.01
```

This can then be run in a bash shell by running
```
bash runs.sh
```
