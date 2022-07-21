#ifndef _INDIVIDUAL_HPP
#define _INDIVIDUAL_HPP

#include <vector>
#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        std::vector < std::vector < double > > hormone_strategy;

        double damage = 0.0;
        int time_since_stressor = 50;

        Individual(int const max_damage
                ,int const max_time_since_stressor
                ,double const hormone_init);

        // copy constructor
        Individual(Individual const &other);

        // birth constructor
        Individual(Individual const &other
                ,Parameters &const params
                ,std::mt19937 &rng_r
                );

        // overloaded assignment operator
        void operator=(Individual const &other);
};

#endif
