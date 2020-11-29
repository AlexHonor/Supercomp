#ifndef _COMPUTATION_CORE_
#define _COMPUTATION_CORE_

#include "Configuration.h"
#include "BlockSolver.h"

class ComputationCore { 
    Configuration conf;
    int rank, world_size;

    public: ComputationCore(INIReader _config) : conf(_config) {}

    public: void Run() {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Vector3Int cell_coord = RankToCell(rank);

        cout << "Process " << rank << "(" << world_size << ") started" << endl 
             << "Coord in grid: { " << cell_coord.x << ", " << cell_coord.y << ", " << cell_coord.z << "}" << endl;

        Vector3Int block_size = conf.cells_num / conf.proc_num;

        Vector3 block_step_size = Vector3::Ones() * (conf.high_border - conf.low_border) * ToVector3(block_size) / ToVector3(conf.cells_num);

        Vector3 bot = ToVector3(cell_coord) * block_step_size;
        Vector3 top = min(ToVector3(cell_coord) * block_step_size, Vector3::Ones() * conf.high_border); 

        Transform transform(block_size * cell_coord, 
                            block_size * cell_coord + block_size,
                            bot,
                            top);

        BlockSolver block(transform, Vector3Int(-1, -1, -1), Vector3Int(-1, -1, -1));

        block.Init();

        //block.Print();
    }

    Vector3Int RankToCell(int rank) {
        return Vector3Int(rank / (conf.proc_num.y * conf.proc_num.z), 
                         (rank % (conf.proc_num.y * conf.proc_num.z)) / conf.proc_num.z,
                         (rank % (conf.proc_num.y * conf.proc_num.z)) % conf.proc_num.z);
    }
};

#endif