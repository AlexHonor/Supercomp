#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>

#include <mpi.h>

#include "external/INIReader.h"

#include "utils.h"
#include "ComputationCore.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc <= 1) {
        cout << "Provide a config file." << endl;

        exit(-1); 
    }
    
    INIReader config(argv[1]);

    if (config.ParseError() != 0) {
        cout << "Unable to open file '" << argv[1] << "'" << endl;
        
        exit(-2); 
    }    

    ComputationCore core(config);

    MPI_Init(&argc, &argv);

    core.Run();

    MPI_Finalize();

    return 0;
}