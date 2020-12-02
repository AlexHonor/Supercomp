#ifndef _BLOCK_SOLVER_
#define _BLOCK_SOLVER_

#include "utils.h"
#include "Matrix.h"
#include "Configuration.h"

#include <map>

double Phi(Vector3 v, Vector3Int size) {
    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z);
}

double Analytical(Vector3 v, Vector3Int size, float t) {
    double alpha_t = M_PI * sqrt(Dot(Vector3(4, 1, 1), ToVector3(size * size)));

    return sin(M_PI / size.x * v.x) * sin(M_PI / size.y * v.y) * sin(M_PI / size.z * v.z) * cos(alpha_t * t + 2 * M_PI);
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

                // @TODO make one pointer
                current = 2;
                prev = 1;
                prevprev = 0;

                n_step = 0;
            }

    public: Matrix& Mat(int step) {
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

//////////////////////////////////////////////////////////////////////////
// Buffers exchange
//////////////////////////////////////////////////////////////////////////

    public: vector<double> AssembleBufferX(int x) {
        vector<double> buffer;
        
        buffer.resize(Mat(0).Size().y * Mat(0).Size().z);

        int i = 0;

        for (int z = 0; z < Mat(0).Size().z; z++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                buffer[i++] = Mat(0)(x, y, z);
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

    public: vector<double> AssembleBufferY(int y) {
        vector<double> buffer;
        
        buffer.resize(Mat(0).Size().x * Mat(0).Size().z);

        int i = 0;

        for (int z = 0; z < Mat(0).Size().z; z++) {
            for (int x = 0; x < Mat(0).Size().x; x++) {
                buffer[i++] = Mat(0)(x, y, z);
            }
        }

        return buffer;
    }

    public: void FillMatrixFromBufferY(vector<double> buffer, int y) {
        int i = 0;
        
        for (int z = 0; z < Mat(0).Size().z; z++) {
            for (int x = 0; x < Mat(0).Size().x; x++) {
                Mat(0)(x, y, z) = buffer[i++];
            }
        }
    }

    public: vector<double> AssembleBufferZ(int z) {
        vector<double> buffer;
        
        buffer.resize(Mat(0).Size().x * Mat(0).Size().y);

        int i = 0;

        for (int y = 0; y < Mat(0).Size().y; y++) {
            for (int x = 0; x < Mat(0).Size().x; x++) {
                buffer[i++] = Mat(0)(x, y, z);
            }
        }

        return buffer;
    }

    public: void FillMatrixFromBufferZ(vector<double> buffer, int z) {
        int i = 0;

        for (int y = 0; y < Mat(0).Size().y; y++) {
            for (int x = 0; x < Mat(0).Size().x; x++) {
                Mat(0)(x, y, z) = buffer[i++];
            }
        }
    }

    public: void InitWithIndices() {
        int i = 0;
        for (int x = 0; x < Mat(-1).Size().x; x++) {
            for (int z = 0; z < Mat(-1).Size().z; z++) {
                for (int y = 0; y < Mat(-1).Size().y; y++) {
                    Mat(-1)(x, y, z) = i++;
                }
            }
        }

        RotateMatrixBuffers();
    } 

    public: void Print() {
        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int z = 0; z < Mat(0).Size().z; z++) {
                for (int y = 0; y < Mat(0).Size().y; y++) {
                    cout << x << "," << y << "," << z << " " << Mat(0)(x, y, z) << endl;
                }
            }
        }
    } 

    static const int BOT_BORDER = 0;
    static const int TOP_BORDER = 1;

    public: void UpdateEdges() {
            vector<MPI_Request> statuses;
            
            vector< vector<double> > incoming_buffers;
            incoming_buffers.resize(6);

            MPI_Request fake;

            int rank = neighbours[hash(Vector3Int(-1, 0, 0))].rank;

            if (rank >= 0) {
                vector<double> outgoing_buffer = AssembleBufferX(0);

                incoming_buffers[0].resize(outgoing_buffer.size());
   
                MPI_Isend(outgoing_buffer.data()    , outgoing_buffer.size()    , MPI_DOUBLE, rank, BOT_BORDER, MPI_COMM_WORLD, &fake);
                MPI_Request_free(&fake);
   
                statuses.push_back(MPI_REQUEST_NULL);

                MPI_Irecv(incoming_buffers[0].data(), incoming_buffers[0].size(), MPI_DOUBLE, rank, TOP_BORDER, MPI_COMM_WORLD, &(*(--statuses.end())));
            }

            rank = neighbours[hash(Vector3Int(1, 0, 0))].rank;

            if (rank >= 0) {
                vector<double> outgoing_buffer = AssembleBufferX(transform.size.x - 1);
   
                incoming_buffers[1].resize(outgoing_buffer.size());
   
                MPI_Isend(outgoing_buffer.data()    , outgoing_buffer.size()    , MPI_DOUBLE, rank, TOP_BORDER, MPI_COMM_WORLD, &fake);
                MPI_Request_free(&fake);

                statuses.push_back(MPI_REQUEST_NULL);

                MPI_Irecv(incoming_buffers[1].data(), incoming_buffers[1].size(), MPI_DOUBLE, rank, BOT_BORDER, MPI_COMM_WORLD, &(*(--statuses.end())));
            }
 
            rank = neighbours[hash(Vector3Int(0, -1, 0))].rank;

            if (rank >= 0) {
                vector<double> outgoing_buffer = AssembleBufferY(0);
   
                incoming_buffers[2].resize(outgoing_buffer.size());
   
                MPI_Isend(outgoing_buffer.data()    , outgoing_buffer.size()    , MPI_DOUBLE, rank, BOT_BORDER, MPI_COMM_WORLD, &fake);
                MPI_Request_free(&fake);

                statuses.push_back(MPI_REQUEST_NULL);

                MPI_Irecv(incoming_buffers[2].data(), incoming_buffers[2].size(), MPI_DOUBLE, rank, TOP_BORDER, MPI_COMM_WORLD, &(*(--statuses.end())));
            }

            rank = neighbours[hash(Vector3Int(0, 1, 0))].rank;

            if (rank >= 0) {
                vector<double> outgoing_buffer = AssembleBufferY(transform.size.y - 1);
   
                incoming_buffers[3].resize(outgoing_buffer.size());
   
                MPI_Isend(outgoing_buffer.data()    , outgoing_buffer.size()    , MPI_DOUBLE, rank, TOP_BORDER, MPI_COMM_WORLD, &fake);
                MPI_Request_free(&fake);

                statuses.push_back(MPI_REQUEST_NULL);

                MPI_Irecv(incoming_buffers[3].data(), incoming_buffers[3].size(), MPI_DOUBLE, rank, BOT_BORDER, MPI_COMM_WORLD,  &(*(--statuses.end())));
            }

            rank = neighbours[hash(Vector3Int(0, 0, -1))].rank;

            if (rank >= 0) {
                vector<double> outgoing_buffer = AssembleBufferZ(0);
   
                incoming_buffers[4].resize(outgoing_buffer.size());
   
                MPI_Isend(outgoing_buffer.data()    , outgoing_buffer.size()    , MPI_DOUBLE, rank, BOT_BORDER, MPI_COMM_WORLD, &fake);
                MPI_Request_free(&fake);

                statuses.push_back(MPI_REQUEST_NULL);

                MPI_Irecv(incoming_buffers[4].data(), incoming_buffers[4].size(), MPI_DOUBLE, rank, TOP_BORDER, MPI_COMM_WORLD,  &(*(--statuses.end())));
            }

            rank = neighbours[hash(Vector3Int(0, 0, 1))].rank;

            if (rank >= 0) {
                vector<double> outgoing_buffer = AssembleBufferZ(transform.size.z - 1);
   
                incoming_buffers[5].resize(outgoing_buffer.size());
   
                MPI_Isend(outgoing_buffer.data()    , outgoing_buffer.size()    , MPI_DOUBLE, rank, TOP_BORDER, MPI_COMM_WORLD, &fake);
                MPI_Request_free(&fake);
   
                statuses.push_back(MPI_REQUEST_NULL);

                MPI_Irecv(incoming_buffers[5].data(), incoming_buffers[5].size(), MPI_DOUBLE, rank, BOT_BORDER, MPI_COMM_WORLD,  &(*(--statuses.end())));
            }

            for (int i = 0; i < statuses.size(); i++) {
                MPI_Wait(&statuses[i], MPI_STATUS_IGNORE);
            }

            if (neighbours[hash(Vector3Int(-1, 0, 0))].rank >= 0) {
                FillMatrixFromBufferX(incoming_buffers[0], -1);
            }

            if (neighbours[hash(Vector3Int(1, 0, 0))].rank >= 0) {
                FillMatrixFromBufferX(incoming_buffers[1], transform.size.x);
            } 

            if (neighbours[hash(Vector3Int(0, -1, 0))].rank >= 0) {
                FillMatrixFromBufferY(incoming_buffers[2], -1);
            } 

            if (neighbours[hash(Vector3Int(0, 1, 0))].rank >= 0) {
                FillMatrixFromBufferY(incoming_buffers[3], transform.size.y);
            }

            if (neighbours[hash(Vector3Int(0, 0, -1))].rank >= 0) {
                FillMatrixFromBufferZ(incoming_buffers[4], -1);
            }

            if (neighbours[hash(Vector3Int(0, 0, 1))].rank >= 0) {
                FillMatrixFromBufferZ(incoming_buffers[5], transform.size.z);
            }            

            //MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            //    int dest, int sendtag,
            //    void *recvbuf, int recvcount, MPI_Datatype recvtype,
            //    int source, int recvtag, MPI_Comm comm, MPI_Status * status)
    }

    public: void Init() {
        #pragma omp parallel
        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                for (int z = 0; z < Mat(0).Size().z; z++) {
                    Vector3 wc = transform.LocalToWorld(Vector3Int(x, y, z));
               
                    Mat(0)(x, y, z) = Phi(wc, Vector3Int::Ones() * conf.cells_num);
                }
            }
        }

        n_step = 0;
        
        ZeroBorders();
        
        GetMaxError();
        
        UpdateEdges();

        RotateMatrixBuffers();

        #pragma omp parallel
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

        ZeroBorders();
        
        n_step++;

        GetMaxError();
    }

    public: double GetMaxError() {
        //cout << "Block " << transform.position << ", " << transform.size << endl; 
        //cout << "Transform" << transform.bot << ", " << transform.top << endl;

        double max_error = 0.0f;

        for (int x = 0; x < Mat(0).Size().x; x++) {
            for (int y = 0; y < Mat(0).Size().y; y++) {
                for (int z = 0; z < Mat(0).Size().z; z++) {
                    double error = fabs(Mat(0)(x, y, z) - Analytical(transform.LocalToWorld(x,y,z), Vector3Int::Ones() * conf.cells_num, n_step * conf.time_step));
                    
                    //cout << Mat(0)(x, y, z) << " " << Analytical(transform.LocalToWorld(x,y,z), transform.size, n_step * conf.time_step) << endl;
                    
                    max_error = max(error, max_error);
                }
            }
        }

        double max_error_total = 666;
        
        MPI_Reduce(&max_error, &max_error_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        if (conf.rank == 0) {
            conf.fout << max_error_total << ",";
        }

        return max_error;
    }

    public: void Iterate() {                
        UpdateEdges();

        #pragma omp parallel
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

        ZeroBorders();
        
        GetMaxError();
    }

    private: void ZeroBorders() {
        /*
        for (int x = 0; x < Mat(0).Size().x; x++) {
            int y = 0;
            
            for (int z = 0; z < Mat(0).Size().z; z++) {
                Mat(0)(x, y, z) = 0;
            }
        }
        
        for (int x = 0; x < Mat(0).Size().x; x++) {
            int z = 0;

            for (int y = 0; y < Mat(0).Size().y; y++) {
                Mat(0)(x, y, z) = 0;    
            }
        }
        */
    } 

    public: void SendResultingMatrix() {
        if (conf.rank != 0) {
            MPI_Send(Mat(0).data.data(), Mat(0).data.size(), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        }
    }
};

#endif
