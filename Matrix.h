#ifndef _MATRIX_
#define _MATRIX_

#include "utils.h"

class Matrix {
    private: Vector3Int size;
    private: vector<double> data;

    public: Matrix(Vector3Int _size) : size(_size + Vector3Int::Ones() * 2) {
        data.resize((_size.x) * (_size.y) * (_size.z));
    }

    public: Vector3Int Size() {
        return size - Vector3Int::Ones() * 2;
    }

    public: bool IsInside(Vector3Int _position) {
        return _position.x >= 0 && _position.y >= 0 && _position.z >= 0 && 
               _position.x < size.x && _position.y < size.y && _position.z < size.z;
    }

    public: double& get(Vector3Int _position) {
        return data[Index(_position)];
    }

    private: size_t Index(Vector3Int _position) {
        return ((_position.x + 1) * size.y + _position.y + 1) * size.z + _position.z + 1;
    }
};

#endif