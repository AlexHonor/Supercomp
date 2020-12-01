#ifndef _CONFIGURATION_
#define _CONFIGURATION_

#include "utils.h"
#include "external/INIReader.h"

class Configuration {
    public: int cells_num;
    public: Vector3Int proc_num;
    public: double time_step, low_border, high_border;

    public: int total_steps;
    public: double total_time;

    public: double grid_step;

    public: int rank;

    public: Configuration(INIReader config) {
        cells_num = config.GetInteger("Computation", "CellsNumber", 32);
        
        total_time  = config.GetReal("Computation", "TotalTime", 1.0);
        total_steps = config.GetInteger("Computation", "IterationNumber", 20);
        
        high_border = config.GetReal("Computation", "HighBorder", 1.0);
        low_border  = config.GetReal("Computation", "LowBorder", 0.0);
    
        grid_step = (high_border - low_border) / cells_num;
    
        time_step   = total_time / total_steps;
    }
};

#endif