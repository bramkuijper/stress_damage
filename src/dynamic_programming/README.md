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
Parameters can also be varied by using a bunch of `for`-loops provided in the Python3 script `vary_params.py`. 
