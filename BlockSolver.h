#ifndef _BLOCK_SOLVER_
#define _BLOCK_SOLVER_

#include "utils.h"
#include "Matrix.h"
#include "Configuration.h"

double Phi(Vector3 v, Vector3Int size) {
    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z);
}

double Analytical(Vector3 v, Vector3Int size, float t) {
    double alpha_t = M_PI * sqrt(Dot(Vector3(1, 1, 1), ToVector3(size * size)));

    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z) * cos(alpha_t * t);
}

struct Neighbour {
    Vector3Int shift;
    int rank;
};

class BlockSolver {
    vector<Matrix> data;
    
    int current, prev, prevprev;

    Transform transform;

    map<int, Neighbour> neighbours;

    Configuration conf;

    public: int n_step;

    public: BlockSolver(Transform _transform, map<int, Neighbour> _neighbours, Configuration _conf) : 
            neighbours(_neighbours),
            transform(_transform),
            conf(_conf) {
                data.push_back(Matrix(_transform.size));
                data.push_back(Matrix(_transform.size));
                data.push_back(Matrix(_transform.size));
                
                current = 2;
                prev = 1;
                prevprev = 0;

                n_step = 0;
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

    public: vector<double> AssembleBufferX(int x) {
        vector<double> buffer;
        
        buffer.resize(Mat(0).Size().y * Mat(0).Size().z);

        for (int z = 0; z < Mat(0).Size().z; z++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                buffer.push_back(Mat(0)(x, y, z));
            }
        }

        return buffer;
    }

    public: void FillMatrixFromBufferX(vector<double> buffer, int x) {
        int i = 0;

        for (int z = 0; z < Mat(0).Size().z; z++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                Mat(0)(x, y, z) = buffer[i++];
            }
        }
    }

    public: void UpdateEdges() {

            {
                AssembleBufferX(0);
                
                int rank = neighbours[hash(Vector3Int(-1, 0, 0))].rank;

               // MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               //              int dest, int sendtag,
               //              void *recvbuf, int recvcount, MPI_Datatype recvtype,
               //              int source, int recvtag, MPI_Comm comm, MPI_Status * status)
            }

            AssembleBufferX(-1);

            //MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            //    int dest, int sendtag,
            //    void *recvbuf, int recvcount, MPI_Datatype recvtype,
            //    int source, int recvtag, MPI_Comm comm, MPI_Status * status)
    }

    public: void Init() {
        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                for (int z = 0; z < Mat(0).Size().z; z++) {
                    Vector3 wc = transform.LocalToWorld(Vector3Int(x, y, z));
               
                    Mat(0)(x, y, z) = Phi(wc, Vector3Int::Ones() * conf.cells_num);
                }
            }
        }
        
        n_step = 0;
    
        cout << "Initialized" << endl;
        cout << GetMaxError() << endl;
        //UpdateEdges();

        RotateMatrixBuffers();

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

                    double accum = -Phi(transform.LocalToWorld(Vector3Int(x, y, z)), Vector3Int::Ones() * conf.cells_num) * 6.0;
                    
                    for (int i = 0; i < sizeof(samples) / sizeof(*samples); i++) {
                        accum += Phi(samples[i], Vector3Int::Ones() * conf.cells_num);
                    }

                    Mat(0)(x, y, z) += Mat(-1)(x, y, z) + 
                                       accum * pow(conf.time_step, 2.0) / 2.0 / pow(conf.grid_step, 2);
                }
            }
        }
        n_step++;

        cout << "Second initialized" << endl;
        cout << GetMaxError() << endl;
    }

    public: double GetMaxError() {
        //cout << "Block " << transform.position << ", " << transform.size << endl; 
        //cout << "Transform" << transform.bot << ", " << transform.top << endl;

        double max_error = 0.0f;

        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                for (int z = 0; z < Mat(0).Size().z; z++) {
                    double error = abs(Mat(0)(x, y, z) - Analytical(transform.LocalToWorld(x,y,z), transform.size, n_step * conf.time_step));
                    
                    //cout << Mat(0)(x, y, z) << " " << Analytical(transform.LocalToWorld(x,y,z), transform.size, n_step * conf.time_step) << endl;
                    
                    max_error = max(error, max_error);
                }
            }
        }

        return max_error;
    }

    public: void Iterate() {
        for (int x = 0; x < Mat(1).Size().x; x++) {
            for (int y = 0; y < Mat(1).Size().y; y++) {
                for (int z = 0; z < Mat(1).Size().z; z++) {
                    Vector3Int samples[] = {
                        Vector3Int(x - 1, y, z),
                        Vector3Int(x, y - 1, z),
                        Vector3Int(x, y, z - 1),
                        Vector3Int(x + 1, y, z),
                        Vector3Int(x, y + 1, z),
                        Vector3Int(x, y, z + 1)
                    };

                    double accum = - Mat(0)(x, y, z) * 4.0 - Mat(-1)(x, y, z);
                    
                    for (int i = 0; i < sizeof(samples) / sizeof(*samples); i++) {
                        accum += Mat(0)(samples[i]);
                    }

                    Mat(1)(x, y, z) = Mat(0)(x, y, z) +
                                        accum * pow(conf.time_step * n_step, 2.0) / pow(conf.grid_step, 2);
                }
            }
        }

        RotateMatrixBuffers();
    
        cout << GetMaxError() << endl;
    }
};

#endif
