#ifndef _COMPUTATION_CORE_
#define _COMPUTATION_CORE_

#include "Configuration.h"
#include "BlockSolver.h"

class ComputationCore { 
    Configuration conf;
    int rank, world_size;

    public: ComputationCore(int argc, char **argv) : conf(argc, argv) {}


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

        //cout << ss.str() << endl;

        return neighbour_ranks;
    }

    Transform CreateBlockTransform(int _rank) {
        Vector3Int cell_coord = RankToCell(_rank);

        Vector3Int block_size = (Vector3Int::Ones() * conf.cells_num + conf.proc_num - Vector3Int::Ones()) / conf.proc_num;

        Vector3 block_step_size = Vector3::Ones() * (conf.high_border - 0) * ToVector3(block_size) / (Vector3::Ones() * conf.cells_num);

        Vector3 bot = ToVector3(cell_coord) * block_step_size;
        Vector3 top = min(ToVector3(cell_coord + Vector3Int::Ones()) * block_step_size, Vector3::Ones() * conf.high_border); 

        return Transform(block_size * cell_coord, 
                         min(block_size * (cell_coord + Vector3Int::Ones()), Vector3Int::Ones() * conf.cells_num) - block_size * cell_coord,
                         bot, top);
    }

    ofstream out;

    public: void Run() {
        double startTime = MPI_Wtime();

        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        // @TODO Refactor
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        conf.rank = rank;

        
        if (rank == 0) {
            //conf.fout << "TotalTime" << ",";
            //conf.fout << "IterationNumber" << ",";
            //conf.fout << "CellsNumber" << ",";
            //conf.fout << "HighBorder" << ",";

            for (int i = 0; i < conf.total_steps; i++) {
            //    conf.fout << i << ",";
            }

            //conf.fout << "TimeSpent" << endl;

            conf.fout << world_size << ",";
            conf.fout << conf.total_time << ",";
            conf.fout << conf.total_steps << ",";
            conf.fout << conf.cells_num << ",";
            conf.fout << conf.high_border << ",";
        }

        conf.proc_num = Vector3Int::Zeros();

        MPI_Dims_create(world_size, 3, (int*)&conf.proc_num);

        //if (conf.rank == 0) {
        //    cout << "Macro grid " << conf.proc_num.x << " - " << conf.proc_num.y << " - " << conf.proc_num.z << endl; 
        //}

        // cout << "Process " << rank << "(" << world_size << ") started" << endl 
        //      << "Coord in grid: { " << cell_coord.x << ", " << cell_coord.y << ", " << cell_coord.z << "}" << endl;

        BlockSolver block(CreateBlockTransform(rank), GetNeighbours(RankToCell(rank)), conf);

        //block.InitWithIndices();

        //block.Print();

        //block.UpdateEdges();

        //block.Print();

        block.Init();

        for (;block.n_step < conf.total_steps; block.n_step++) {
            block.Iterate();
        }

        block.SendResultingMatrix();

        if (rank == 0) {
            AssembleResult(block);

            conf.fout << MPI_Wtime() - startTime << endl;
        }
    }

    Matrix AssembleResult(BlockSolver block) {
        Matrix result(Vector3Int::Ones() * conf.cells_num);

        vector<MPI_Request> statuses;
        statuses.resize(world_size);

        vector<vector<double> > blocks;
        blocks.resize(world_size);
        
        for (int i = 1; i < world_size; i++) {
            Transform block_transform = CreateBlockTransform(i);

            statuses.push_back(MPI_REQUEST_NULL);
            blocks[i].resize((block_transform.size.x + 2) * (block_transform.size.y + 2) * (block_transform.size.z + 2));
         
            MPI_Irecv(blocks[i].data(), blocks[i].size(), MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &statuses[i]);
        }
        
        for (int i = 1; i < statuses.size(); i++) {
            MPI_Wait(&statuses[i], MPI_STATUS_IGNORE);
        }
        
        for (int i = 1; i < world_size; i++) {
            Transform block_transform = CreateBlockTransform(i);

            Matrix local_matrix(block_transform.size);
            
            for (int j = 0; j < blocks[i].size(); j++) {
                local_matrix.data[j] = blocks[i][j];
            }
            for (int x = 0; x < local_matrix.Size().x; x++) {
                for (int y = 0; y < local_matrix.Size().y; y++) {
                    for (int z = 0; z < local_matrix.Size().z; z++) {
                        result(block_transform.LocalToWorldIndex(Vector3Int(x, y, z))) = local_matrix(x, y, z);
                    }
                }
            }
        }

        for (int x = 0; x < block.Mat(0).Size().x; x++) {
            for (int y = 0; y < block.Mat(0).Size().y; y++) {
                for (int z = 0; z < block.Mat(0).Size().z; z++) {
                    result(Vector3Int(x, y, z)) = block.Mat(0)(x, y, z);
                }
            }
        }
        
        stringstream ss;
        ss << conf.cells_num << "-" << conf.time_step << ".bin";

        dump_vector(result.data, ss.str().c_str());

        return result;
    }

    static void dump_vector(const std::vector<double> &v, const char* filename) {
        std::ofstream out_file;
        out_file.open(filename, std::ios::out | std::ios::binary);
        if (!out_file) {
            std::cerr << "can't open file " << filename << std::endl;
            return;
        }
        out_file.write((char *) v.data(), v.size()*sizeof(double));
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