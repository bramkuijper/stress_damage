#include <vector>
#include <random>
#include <algorithm>
#include "parameters.hpp"
#include "individual.hpp"

Individual::Individual(int const max_damage, 
        int const max_time_since_stressor,
        double const hormone_init
        ) :
    damage{0.0}
    ,time_since_stressor{max_time_since_stressor}
    ,hormone_strategy(max_damage, std::vector < double >(max_time_since_stressor + 1, hormone_init))
{
} // end Individual()

// copy constructor - make 1 individual with another
Individual::Individual(Individual const &other) :
    damage{other.damage}
    ,time_since_stressor{other.max_time_since_stressor}
    ,hormone_strategy{other.hormone_strategy}
{
} // end copy constructor

// birth constructor - make 1 individual with another
Individual::Individual(
        Individual const &other
        ,Parameters &const params
        ,std::mt19937 &rng_r
        ) :
    damage{other.damage}
    ,time_since_stressor{other.max_time_since_stressor}
    ,hormone_strategy{other.hormone_strategy}
{
    // setup distributions for the mutation
    std::uniform_real_distribution <double> uniform(0,1.0);
    std::normal_distribution <double> normal(0,params.sdmu);

    double h;

    // go through all the strategy loci and mutate them
    for (int damage_idx = 0; damage_idx < max_damage; ++damage_idx)
    {
        for (int time_since_stressor_idx = 0; 
                time_since_stressor_idx <= max_time_since_stressor; 
                ++time_since_stressor_idx)
        {
            if (uniform(rng_r) < params.mu_hormone)
            {
                h = hormone_strategy[damage_idx][time_since_stressor_idx] + 
                    normal(rng_r);

                hormone_strategy[damage_idx][time_since_stressor_idx] = 
                    std::clamp(h, 0.0, 1.0);
            }
        }// end for time_since_stressor_idx
    } // end for int damage_idx
} // end birth constructor




// assignment operator make 1 individual with another
void Individual::operator=(Individual const &other)
{
    damage = other.damage;
    time_since_stressor = other.time_since_stressor;
    hormone_strategy = other.hormone_strategy;
} // end copy constructor
