#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>

#include <mpi.h>

#include "external/INIReader.h"

using namespace std;

struct Vector3 {
    double x, y, z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

struct Vector3Int {
    int32_t x, y, z;

    Vector3Int() : x(0), y(0), z(0) {}
    Vector3Int(int32_t _x, int32_t _y, int32_t _z) : x(_x), y(_y), z(_z) {}
};

ostream& operator<<(ostream &os, const Vector3Int &v)
{
    os << "{ " << v.x << ", " << v.y << ", " << v.z << "}";
    return os;
}

ostream& operator<<(ostream &os, const Vector3 &v)
{
    os << "{ " << v.x << ", " << v.y << ", " << v.z << "}";
    return os;
}

double FunctionPhi(Vector3 v, Vector3Int size) {
    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z);
}

struct Transform {
    public: Vector3Int position, size;
    public: Vector3 bot, top;

    public: Transform(Vector3Int p, Vector3Int s, Vector3 b, Vector3 t) : position(p), size(s), bot(b), top(t) {}

    public: Vector3 LocalToWorld(Vector3Int local_position) {
        Vector3 result;
        
        result.x = (position.x + local_position.x) / (double)size.x * (top.x - bot.x);
        result.y = (position.y + local_position.y) / (double)size.y * (top.y - bot.y);
        result.z = (position.z + local_position.z) / (double)size.z * (top.z - bot.z);
        
        return result;
    }
};

class Matrix {
    private: Vector3Int size;
    private: vector<double> data;

    public: Matrix(Vector3Int _size) : size(_size) {
        data.resize(_size.x * _size.y * _size.z);
    }

    public: Vector3Int Size() {
        return size;
    }

    public: bool IsInside(Vector3Int _position) {
        return _position.x >= 0 && _position.y >= 0 && _position.z >= 0 && 
               _position.x < size.x && _position.y < size.y && _position.z < size.z;
    }

    public: double& get(Vector3Int _position) {
        return data[Index(_position)];
    }

    private: size_t Index(Vector3Int _position) {
        return (_position.x * size.y + _position.y) * size.z + _position.z;
    }
};

class BlockSolver {
    Matrix matrix;

    Transform transform;

    Vector3Int low_neighbours, high_neighbours;

    public: BlockSolver(Transform _transform, Vector3Int _low_neighbours, Vector3Int _high_neighbours ) : matrix(_transform.size), low_neighbours(_low_neighbours), high_neighbours(_high_neighbours), transform(_transform) {}

    public: void Init() {
        for (int x = 0; x < matrix.Size().x; x++) {
            for (int y = 0; y < matrix.Size().y; y++) {
                for (int z = 0; z < matrix.Size().z; z++) {
                    Vector3 world_coord = transform.LocalToWorld(Vector3Int(x, y, z));

                    matrix.get(Vector3Int(x, y, z)) = FunctionPhi(world_coord, transform.size);
                }
            }
        }
    }

    public: void Print() {
        cout << "Block " << transform.position << ", " << transform.position << endl; 
        cout << "Transform" << transform.bot << ", " << transform.bot << endl;
                
        for (int x = 0; x < matrix.Size().x; x++) {
            for (int y = 0; y < matrix.Size().y; y++) {
                for (int z = 0; z < matrix.Size().z; z++) {
                    cout << matrix.get(Vector3Int(x, y, z)) << " ";
                }
            }
        }
    }

    public: void Eval() {

    }
};

class ComputationCore { 
    INIReader config;
    int rank, world_size;

    public: ComputationCore(INIReader _config) : config(_config) {}

    public: void Run() {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        cout << "Process " << rank << "(" << world_size << ") started" << endl;        
    
        int x_cells = 1;
        int y_cells = 3;
        int z_cells = 2;

        double low_border = 0.0;
        double high_border = 1.0;

        int grid_size = 4;

        double grid_step_size = (high_border - low_border) / grid_size;
        
        int block_size = grid_size / world_size;

        int bot = rank * block_size;
        int top = min(bot + block_size, grid_size);

        Transform transform(Vector3Int(bot, 0, 0), 
                            Vector3Int(top, grid_size, grid_size),
                            Vector3(bot * grid_step_size, low_border, low_border),
                            Vector3(top * grid_step_size, high_border, high_border));

        BlockSolver block(transform, Vector3Int(-1, -1, -1), Vector3Int(-1, -1, -1));

        block.Init();

        block.Print();
    }
};

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