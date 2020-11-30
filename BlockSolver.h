#ifndef _BLOCK_SOLVER_
#define _BLOCK_SOLVER_

#include "utils.h"
#include "Matrix.h"
#include "Configuration.h"

double Phi(Vector3 v, Vector3Int size) {
    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z);
}

class BlockSolver {
    Matrix a, b, c;
    vector<Matrix> data;
    
    int current, prev, prevprev;

    Transform transform;

    Vector3Int low_neighbours, high_neighbours;

    Configuration conf;

    public: BlockSolver(Transform _transform, Vector3Int _low_neighbours, Vector3Int _high_neighbours, Configuration _conf) : 
            low_neighbours(_low_neighbours), 
            high_neighbours(_high_neighbours), 
            transform(_transform),
            conf(_conf),
            a(transform.size),
            b(transform.size),
            c(transform.size) {
                data.push_back(a);
                data.push_back(b);
                data.push_back(c);
                
                current = 2;
                prev = 1;
                prevprev = 0;
            }

    private: Matrix& Mat(int step) {
        int index[] = {
            current,
            prev,
            prevprev
        };
        return data[index[step + 1]];
    } 

    private: void RotateMatrixBuffers() {
        prev = (prev + 1) % 3;
        current = (current + 1) % 3;
        prevprev = (prevprev + 1) % 3;
    } 

    public: void Init() {
        for (int x = 0; x < Mat(-1).Size().x; x++) {
            for (int y = 0; y < Mat(-1).Size().y; y++) {
                for (int z = 0; z < Mat(-1).Size().z; z++) {
                    Vector3 wc = transform.LocalToWorld(Vector3Int(x, y, z));
                    
                    //cout << wc << " ";
                
                    Mat(-1).get(Vector3Int(x, y, z)) = Phi(wc, transform.size);
                }
            }
        }

        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                for (int z = 0; z < Mat(0).Size().z; z++) {
                    Vector3 samples[] = {
                        transform.LocalToWorld(Vector3Int(x - 1, y, z)),
                        transform.LocalToWorld(Vector3Int(x, y - 1, z)),
                        transform.LocalToWorld(Vector3Int(x, y, z - 1)),
                        transform.LocalToWorld(Vector3Int(x + 1, y, z)),
                        transform.LocalToWorld(Vector3Int(x, y + 1, z)),
                        transform.LocalToWorld(Vector3Int(x, y, z + 1))
                    };

                    double accum = Phi(transform.LocalToWorld(Vector3Int(x, y, z)), transform.size) * 6.0;
                    
                    for (int i = 0; i < sizeof(samples) / sizeof(*samples); i++) {
                        accum += Phi(samples[i], transform.size);
                    }

                    Mat(0).get(Vector3Int(x, y, z)) += Mat(-1).get(Vector3Int(x, y, z)) + 
                                                       Phi(transform.LocalToWorld(Vector3Int(x, y, z)), transform.size) +
                                                       accum * pow(conf.time_step, 2.0) / 2.0 / pow(conf.grid_step, 2);
                }
            }
        }
    }

    public: void Print() {
        cout << "Block " << transform.position << ", " << transform.size << endl; 
        cout << "Transform" << transform.bot << ", " << transform.top << endl;
                
        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                for (int z = 0; z < Mat(0).Size().z; z++) {
                    cout << Mat(0).get(Vector3Int(x, y, z)) << " ";
                }
            }
        }
    }

    public: void Eval() {
        for (int x = 0; x < Mat(1).Size().x; x++) {
            for (int y = 0; y < Mat(1).Size().y; y++) {
                for (int z = 0; z < Mat(1).Size().z; z++) {
                    Vector3 samples[] = {
                        transform.LocalToWorld(Vector3Int(x - 1, y, z)),
                        transform.LocalToWorld(Vector3Int(x, y - 1, z)),
                        transform.LocalToWorld(Vector3Int(x, y, z - 1)),
                        transform.LocalToWorld(Vector3Int(x + 1, y, z)),
                        transform.LocalToWorld(Vector3Int(x, y + 1, z)),
                        transform.LocalToWorld(Vector3Int(x, y, z + 1))
                    };

                    double accum = Phi(transform.LocalToWorld(Vector3Int(x, y, z)), transform.size) * 6.0;
                    
                    for (int i = 0; i < sizeof(samples) / sizeof(*samples); i++) {
                        accum += Phi(samples[i], transform.size);
                    }

                    Mat(1).get(Vector3Int(x, y, z)) += Phi(transform.LocalToWorld(Vector3Int(x, y, z)), transform.size) +
                                                       accum * pow(conf.time_step, 2.0) / 2.0 / pow(conf.grid_step, 2);
                }
            }
        }
    }
};

#endif
