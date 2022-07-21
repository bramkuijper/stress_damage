#ifndef _INDIVIDUAL_HPP
#define _INDIVIDUAL_HPP

#include <vector>
#include "parameters.hpp"

class Individual
{
    public:
        std::vector < std::vector < double > > hormone_strategy;

        double damage = 0.0;

        Individual();

        // copy constructor
        Individual(Individual const &other);

        // overloaded assignment operator
        void operator=(Individual const &other);
};

#endif
