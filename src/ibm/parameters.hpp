#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

class Parameters
{
    public:
        double Kmort = 0.001;
        double Kfec = 0.001;
        
        long unsigned max_time = 100;

        // duration of the season
        unsigned season_duration = 10;

        // onset of breeding season
        // where breeding_t_start < season_duration
        unsigned breeding_t_start = 5;

        // end of breeding season
        unsigned breeding_t_end = 5;
};

#endif
