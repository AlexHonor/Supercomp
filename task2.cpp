#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>

#include <mpi.h>

#include "utils.h"
#include "ComputationCore.h"

ofstream Configuration::fout;

using namespace std;

int main(int argc, char** argv) {
    ComputationCore core(argc, argv);

    MPI_Init(&argc, &argv);

    core.Run();

    MPI_Finalize();

    return 0;
}