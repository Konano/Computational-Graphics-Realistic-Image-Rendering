/*
* @Author: Konano
* @Date:   2019-06-21 13:31:34
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:31:34
*/

#ifndef __VEC3_H__
#define __VEC3_H__

#include "utils.hpp"

struct Vec {
    double x, y, z;
    Vec(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z) {}
    Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec &operator+=(const Vec &b) { x += b.x, y += b.y, z += b.z; return *this;}
    Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec &operator-=(const Vec &b) { x -= b.x, y -= b.y, z -= b.z; return *this;}
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec operator/(double b) const { return Vec(x / b, y / b, z / b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec cross(const Vec &b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    Vec& norm() { return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
    Vec& clp() { x=clamp(x); y=clamp(y), z=clamp(z); return *this; }
    double len() const { return sqrt(x*x + y*y + z*z); }
    double len2() const { return x*x + y*y + z*z; }
    double max() const { return (x>y && x>z) ? x : (y>z ? y : z); }
    double min() const { return (x<y && x<z) ? x : (y<z ? y : z); }
    void print() const { printf("%.12lf %.12lf %.12lf\n", x, y, z); }
    bool smaller(const Vec& b) const { return x<=b.x+eps && y<=b.y+eps && z<=b.z+eps; }
    bool bigger(const Vec& b) const { return x>=b.x-eps && y>=b.y-eps && z>=b.z-eps; }
};

Vec Vmax(const Vec &a, const Vec &b) {
    return Vec(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

Vec Vmin(const Vec &a, const Vec &b) {
    return Vec(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

#endif // __VEC3_H__