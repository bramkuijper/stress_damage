#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <vector>
#include "parameters.hpp"
#include "individual.hpp"

class Simulation
{
    public:
        Parameters params;
        std::vector <Individual> Pop;

        Simulation(Parameters const &params);
   
        // actually run the simulation
        void run();
};

#endif
