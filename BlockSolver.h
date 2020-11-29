#ifndef _BLOCK_SOLVER_
#define _BLOCK_SOLVER_

#include "utils.h"
#include "Matrix.h"

double FunctionPhi(Vector3 v, Vector3Int size) {
    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z);
}

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

#endif
