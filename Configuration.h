#ifndef _CONFIGURATION_
#define _CONFIGURATION_

#include "utils.h"
#include "external/INIReader.h"

class Configuration {
    public: Vector3Int cells_num;
    public: Vector3Int proc_num;
    public: double time_step, low_border, high_border;

    public: Configuration(INIReader config) {
        cells_num.x = config.GetInteger("Computation", "CellsNumberX", 32);
        cells_num.y = config.GetInteger("Computation", "CellsNumberY", 32);
        cells_num.z = config.GetInteger("Computation", "CellsNumberZ", 32);
        
        proc_num.x = config.GetInteger("Computation", "ProcNumX", 1);
        proc_num.y = config.GetInteger("Computation", "ProcNumY", 1);
        proc_num.z = config.GetInteger("Computation", "ProcNumZ", 1);
        
        time_step   = config.GetReal("Computation", "TimeStep", 0.1);
        high_border = config.GetReal("Computation", "HighBorder", 1.0);
        low_border  = config.GetReal("Computation", "LowBorder", 0.0);
    }
};

#endif