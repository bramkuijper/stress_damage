#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#include <string>

class Parameters
{
    public:
        double Kmort = 0.001;
        double Kfec = 0.001;

        double pLeave = 0.5;
        double pArrive = 0.5;
        
        long unsigned max_time = 100;

        // onset of breeding season
        // where breeding_t_start < season_duration
        unsigned breeding_t_start = 5;

        // end of breeding season
        unsigned breeding_t_end = 5;
        int max_damage = 50;
        int max_time_since_stressor = 10;
        double hormone_init = 0.3;

        unsigned int n = 1000;

        double baseline_mortality = 0.002;
        double pkill = 0.5;

        std::string base_name; 
};

#endif
