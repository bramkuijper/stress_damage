#include <vector>
#include <random>
#include <string>
#include "parameters.hpp"
#include "individual.hpp"
#include "simulation.hpp"

Simulation::Simulation(Parameters const &params) :
    params{params}
    ,rd{} // initialize random device, see *.h file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,data_file{params.base_name.c_str()} // initialize the data file by giving it a name
    ,Pop{params.n, Individual(params.max_damage
            ,params.max_time_since_stressor
            ,params.hormone_init)}
{}

void Simulation::run()
{
    // each individual can expereince
    for (int time_step = 0; 
            time_step <= params.max_time; ++time_step)
    {
        survive();
        reproduce();

        if (time_step % data_interval == 0)
        {
            write_data();
        }
    }

    write_parameters();
} // end Simulation::run()

void Simulation::survive()
{
    for (int individual_idx = 0; individual_idx < params.n; ++individual_idx)
    {
        double survival_probability = 1.0;

        // if currently experiencing stressor
        if (Pop[individual_idx].time <= 1)
        {
            if (uniform(rng_r) < params.pLeave)
            {
                survival_probability *= 1.0 - params.background_mortality;


    }
} // end Simulation::survive()
