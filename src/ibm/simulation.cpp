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
    for (time_step = 0; 
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

void Simulation::reproduce()
{
    if (Pop.size() < 1)
    {
        write_parameters();
        exit(1);
    }

    int n_recruits = params.n - Pop.size();

    assert(n_recruits >= 0);

    // sample parents and reproduce
    unsigned int random_parent;

    std::uniform_int_distribution < int > parent_sampler(0, Pop.size() - 1);

    for (int sample_idx = 0; sample_idx < n_recruits; ++sample_idx)
    {
        random_parent = parent_sampler(rng_r);

        // add kid (using the birth constructor in individual.cpp)
        Pop.push_back(
                Individual(
                    Pop[random_parent]
                    ,params
                    ,rng_r)
                );
    }// end int sample_idx

} // end void reproduce()


void Simulation::survive()
{
    // auxiliary variables
    double prob_survive, hormone_i, prob_stressor_present; 
    int damage_i;
    int time_since_stressor_i;

    // go through all individuals and assess survival
    for (int individual_idx = 0; 
            individual_idx < params.n; ++individual_idx)
    {
        // express damage and get time since individual last encountered 
        // stressor
        damage_i = Pop[individual_idx].damage;
        time_since_stressor_i = Pop[individual_idx].time_since_stressor;

        // prob_survive baseline is 1.0 - mu
        prob_survive = 1.0 - std::min(1.0, params.background_mortality + 
            params.Kmort * damage_i);

        // calculate prob that stressor is currently present
        prob_stressor_present = time_since_stressor_i <= 1 ?
            1.0 - params.pLeave // stressor still around, will it stay?
            :
            params.pArrive; // stressor not around, will it appear?
   
        // if stressor present and attacks
        if (uniform(rng_r) < prob_stressor_present && 
                uniform(rng_r) < params.p_attack)
        {
            // express damage
            hormone_i = Pop[individual_idx].hormone_strategy[damage_i][time_since_stressor_i];

            prob_survive *= 1.0 - prob_killed(hormone_i);
            
            Pop[individual_idx].time_since_stressor = 1;
        }
        else
        {
            ++Pop[individual_idx].time_since_stressor;
        }

        if (uniform(rng_r) < 1.0 - prob_survive)
        {
            Pop.erase(Pop.begin() + individual_idx);
            --individual_idx;
        }
    } // end for (int individual_idx = 0; 

} // end Simulation::survive()

void Simulation::write_data()
{
    // make vector mean_hormone_strategy[damage_i][time_since_attack_j]
    std::vector < std::vector < double > > mean_hormone_strategy(
            max_damage, std::vector < double >( max_time_since_stressor + 1, 0.0));

    // sum of squares also, to calculate variances
    std::vector < std::vector < double > > ss_hormone_strategy(
            max_damage, std::vector < double >( max_time_since_stressor + 1, 0.0));

    for (unsigned int individual_idx = 0; individual_idx < Pop.size(); ++individual_idx)
    {
        // go through all the strategy loci and mutate them
        for (int damage_idx = 0; damage_idx < max_damage; ++damage_idx)
        {
            for (int time_since_stressor_idx = 0; 
                    time_since_stressor_idx <= max_time_since_stressor; 
                    ++time_since_stressor_idx)
            {
                h = Pop[individual_idx].hormone_strategy[damage_idx][time_since_stressor_idx];
                
                mean_hormone_strategy[damage_idx][time_since_stressor_idx] += h;
                ss_hormone_strategy[damage_idx][time_since_stressor_idx] += h * h;
            }
        }
    } // end invididual_idx

    double mean, var;
    // go through all the strategy loci and mutate them
    for (int damage_idx = 0; damage_idx < max_damage; ++damage_idx)
    {
        for (int time_since_stressor_idx = 0; 
                time_since_stressor_idx <= max_time_since_stressor; 
                ++time_since_stressor_idx)
        {
            mean = mean_hormone_strategy[damage_idx][time_since_stressor_idx] / Pop.size();
            var = ss_hormone_strategy[damage_idx][time_since_stressor_idx] / Pop.size()
                - mean * mean;

            data_file << time_step << ";" 
                << damage_idx << ";" 
                << time_since_stressor_idx << ";" 
                << mean << ";"
                << var << std::endl;
        }
    }
} // end write_data()

void Simulation::write_parameters()
{
    data_file 
        << "Kmort;" << Kmort << std::endl
        << "Kfec;" << Kfec << std::endl
        << "pLeave;" << pLeave << std::endl
        << "pArrive;" << pArrive << std::endl
        << "risk;" << (pArrive / (pArrive + pLeave)) << std::endl
        << "autocorr;" << 1.0 - pArrive - pLeave << std::endl
        << "alpha;" << alpha << std::endl
        << "hormone_init;" << hormone_init << std::endl
        << "baseline_mortality;" << baseline_mortality << std::endl

} // end write_parameters()

void Simulation::write_data_headers()
{
    data_file << "time;damage;time_since_stressor;meanh;varh" << std::endl;

} // end write_parameters()

double prob_killed(double const hormone_level)
{
    return(std::max(0.0, 1.0 - pow(hormone_level, params.alpha)));
} // end pkilled
