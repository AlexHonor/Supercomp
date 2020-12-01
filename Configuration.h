#ifndef _CONFIGURATION_
#define _CONFIGURATION_

#include "utils.h"
#include "external/INIReader.h"

class Configuration {
    public: int cells_num;
    public: Vector3Int proc_num;
    public: double time_step, low_border, high_border;

    public: double grid_step;

    public: Configuration(INIReader config) {
        cells_num = config.GetInteger("Computation", "CellsNumber", 32);
        
        proc_num.x = config.GetInteger("Computation", "ProcNumX", 1);
        proc_num.y = config.GetInteger("Computation", "ProcNumY", 1);
        proc_num.z = config.GetInteger("Computation", "ProcNumZ", 1);
        
        high_border = config.GetReal("Computation", "HighBorder", 1.0);
        low_border  = config.GetReal("Computation", "LowBorder", 0.0);
    
        grid_step = (high_border - low_border) / cells_num;
    
        time_step   = grid_step * grid_step / 2.0;
    }
};

#endif