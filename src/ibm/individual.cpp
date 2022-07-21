#include <vector>
#include "parameters.hpp"
#include "individual.hpp"

Individual::Individual(int const max_damage, 
        int const max_time_since_stressor,
        double const hormone_init
        ) :
    damage{0.0}
    ,time_since_stressor{max_time_since_stressor}
{
    for (int damage_idx = 0; damage_idx < max_damage; ++damage_idx)
    {
        //initialize an array of size max_time_since_stressor
        //with a hormone value of hormone_init
        std::vector < double > sub_array(max_time_since_stressor + 1, hormone_init);

        hormone_strategy.push_back(sub_array);
    }
} // end Individual()

// copy constructor - make 1 individual with another
Individual::Individual(Individual const &other) :
    damage{other.damage}
    ,time_since_stressor{max_time_since_stressor}
    ,hormone_strategy{other.hormone_strategy}
{
} // end copy constructor

// assignment operator make 1 individual with another
void Individual::operator=(Individual const &other)
{
    damage = other.damage;
    time_since_stressor = time_since_stressor.damage;
    hormone_strategy = other.hormone_strategy;
} // end copy constructor