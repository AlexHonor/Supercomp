#ifndef _CONFIGURATION_
#define _CONFIGURATION_

#include <fstream>
#include <cstdlib> 

#include "utils.h"

class Configuration {
    public: int cells_num;
    public: Vector3Int proc_num;
    public: double time_step, high_border;

    public: int total_steps;
    public: double total_time;

    public: double grid_step;

    public: int rank;

    public: static ofstream fout;

    public: Configuration(int argc, char **argv) {
        if (argc < 6) {
            cout << "Provide arguments [CellsNumber] [TotalTime] [IterationNumber] [HighBorder] [OutCSVFilename]" << endl;

            exit(-1); 
        }

        cells_num = atoi(argv[1]);
        total_time  = atof(argv[2]);
        total_steps = atoi(argv[3]);
        high_border = atof(argv[4]);
        
        grid_step = (high_border - 0) / cells_num;
    
        time_step = total_time / total_steps;
    
        fout.open(argv[5], ofstream::out | ofstream::app);
    }
};

#endif