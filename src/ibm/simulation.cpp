#include <vector>
#include "parameters.hpp"
#include "individual.hpp"

Simulation::Simulation(Parameters const &params) :
    params{params}
    ,Pop{params.n, Individual()}
{}

void Simulation::run()
{}
