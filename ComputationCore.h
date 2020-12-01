#ifndef _COMPUTATION_CORE_
#define _COMPUTATION_CORE_

#include "Configuration.h"
#include "BlockSolver.h"

class ComputationCore { 
    Configuration conf;
    int rank, world_size;

    public: ComputationCore(INIReader _config) : conf(_config) {}


    // @TODO move to block solver?
    // wrong neighbours
    // how to perform sending and recieving
    private: map<int, Neighbour> GetNeighbours(Vector3Int block_position) {
        map<int, Neighbour> neighbour_ranks;

        stringstream ss;
        
        vector<Vector3Int> shifts;
        shifts.push_back(Vector3Int(-1, 0, 0));
        shifts.push_back(Vector3Int(+1, 0, 0));
        shifts.push_back(Vector3Int(0, -1, 0));
        shifts.push_back(Vector3Int(0, +1, 0));
        shifts.push_back(Vector3Int(0, 0, -1));
        shifts.push_back(Vector3Int(0, 0, +1));

        ss << "Cell " << block_position << endl;

        for (int i = 0; i < shifts.size(); i++) {
            Vector3Int neighbour_position = block_position + shifts[i];
            int neighbour = CellToRank(neighbour_position);
            
            if (neighbour_position.x == -1) {
                neighbour = CellToRank(Vector3Int(conf.proc_num.x - 1, neighbour_position.y, neighbour_position.z));
            } else if (neighbour_position.x == conf.proc_num.x) {
                neighbour = CellToRank(Vector3Int(0, neighbour_position.y, neighbour_position.z));
            } else if (neighbour_position.y == -1 || neighbour_position.y == conf.proc_num.y ||
                neighbour_position.z == -1 || neighbour_position.z == conf.proc_num.z) {
                neighbour = -1;
            }

            neighbour_ranks[hash(shifts[i])].shift = shifts[i];
            neighbour_ranks[hash(shifts[i])].rank = neighbour;
        }

        return neighbour_ranks;
    }
 

    public: void Run() {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        conf.proc_num = Vector3Int::Zeros();

        MPI_Dims_create(world_size, 3, (int*)&conf.proc_num);

        Vector3Int cell_coord = RankToCell(rank);

        cout << "Process " << rank << "(" << world_size << ") started" << endl 
             << "Coord in grid: { " << cell_coord.x << ", " << cell_coord.y << ", " << cell_coord.z << "}" << endl;

        Vector3Int block_size = Vector3Int::Ones() * conf.cells_num / conf.proc_num;

        Vector3 block_step_size = Vector3::Ones() * (conf.high_border - conf.low_border) * ToVector3(block_size) / (Vector3::Ones() * conf.cells_num);

        Vector3 bot = ToVector3(cell_coord) * block_step_size;
        Vector3 top = min(ToVector3(cell_coord + Vector3Int::Ones()) * block_step_size, Vector3::Ones() * conf.high_border); 

        Transform transform(block_size * cell_coord, 
                            min(block_size * (cell_coord + Vector3Int::Ones()), Vector3Int::Ones() * conf.cells_num) - block_size * cell_coord,
                            bot,
                            top);

        BlockSolver block(transform, GetNeighbours(cell_coord), conf);


        //block.InitWithIndices();

        //block.Print();

        //block.UpdateEdges();

        //block.Print();

        block.Init();

        for (;block.n_step < 20; block.n_step++) {
            block.Iterate();
        }
    }

    Vector3Int RankToCell(int rank) {
        return Vector3Int(rank / (conf.proc_num.y * conf.proc_num.z), 
                         (rank % (conf.proc_num.y * conf.proc_num.z)) / conf.proc_num.z,
                         (rank % (conf.proc_num.y * conf.proc_num.z)) % conf.proc_num.z);
    }

    int CellToRank(Vector3Int cell) {
        return cell.x * (conf.proc_num.y * conf.proc_num.z) +  
               cell.y * conf.proc_num.z + 
               cell.z;
    }
};

#endif