#include "parameters.hpp"
#include "simulation.hpp"

int main(int argc, char **argv)
{
    Parameters parms;

    Simulation sim(parms);

    sim.run();
} // end main
