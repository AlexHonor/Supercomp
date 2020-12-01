#ifndef _UTILS_
#define _UTILS_

#include <iostream>

using namespace std;

template<typename T>
struct TVector3 {
    T x, y, z;

    TVector3() : x(0), y(0), z(0) {}
    TVector3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}

    TVector3 operator+(const TVector3 r) {
        TVector3 result;

        result.x = x + r.x;
        result.y = y + r.y;
        result.z = z + r.z;

        return result;
    }
    
    TVector3 operator-(const TVector3 r) {
        TVector3 result;

        result.x = x - r.x;
        result.y = y - r.y;
        result.z = z - r.z;

        return result;
    }

    TVector3 operator*(const TVector3 r) {
        TVector3 result;

        result.x = x * r.x;
        result.y = y * r.y;
        result.z = z * r.z;

        return result;
    }

    TVector3 operator*(const T r) {
        TVector3 result;

        result.x = x * r;
        result.y = y * r;
        result.z = z * r;

        return result;
    }

    TVector3 operator/(const TVector3 r) {
        TVector3 result;

        result.x = x / r.x;
        result.y = y / r.y;
        result.z = z / r.z;

        return result;
    }

    static TVector3 Ones() {
        return TVector3(1, 1, 1);
    }
    
    static TVector3 Zeros() {
        return TVector3();
    }
};

typedef TVector3<int32_t> Vector3Int;
typedef TVector3<double> Vector3;

template<typename T>
T Dot(TVector3<T> a, TVector3<T> b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 ToVector3(Vector3Int v) {
    return Vector3(v.x, v.y, v.z);
}

Vector3Int ToVector3Int(Vector3 v) {
    return Vector3Int(v.x, v.y, v.z);
}

template<typename T>
TVector3<T> min(TVector3<T> a, TVector3<T> b) {
    return TVector3<T>(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

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

int hash(int x, int y, int z) {
    return (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
}

int hash(Vector3Int v) {
    return hash(v.x, v.y, v.z);
}

struct Transform {
    public: Vector3Int position, size;
    public: Vector3 bot, top;

    public: Transform(Vector3Int p, Vector3Int s, Vector3 b, Vector3 t) : position(p), size(s), bot(b), top(t) {}

    public: Vector3 LocalToWorld(int x, int y, int z) {
        return LocalToWorld(Vector3Int(x, y, z));
    }

    public: Vector3 LocalToWorld(Vector3Int local_position) {
        return ToVector3(position + local_position) / ToVector3(size) * (top - bot);
    }
};

#endif