/*
* @Author: Konano
* @Date:   2019-06-19 15:57:28
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-19 18:38:58
*/

#include <cmath>
#include <cstdlib>
#include <cstdio>

#define PI 3.14159265358979323846

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1.0 / 2.2) * 255 + 0.5); } // Gamma Correction

struct Vec {
    double x, y, z;
    Vec(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z) {}
    Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec cross(const Vec &b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    // Vec& norm() { return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
};

struct Ray {
    Vec ori, dir;
    Ray(Vec o, Vec d) : ori(o), dir(d) {}
};

struct Sphere
{
    double rad;
    Vec pos;
    // Vec position, emission, color;

    Sphere(double r, Vec p) : rad(r), pos(p) {}

    double intersect(const Ray &r) const {
        const double eps = 1e-4;
        Vec dis = r.ori - pos;
        double a = r.dir * r.dir;
        double b = r.dir * dis * 2;
        double c = dis * dis - rad * rad;
        double detal =  b * b - 4 * a * c;
        if (detal < 0) return 0;
        detal = sqrt(detal);
        double t1 = (- b - detal) / 2 / a;
        double t2 = (- b + detal) / 2 / a;
        if (t1 > eps) return t1;
        else if (t2 > eps) return t2;
        else return 0;
    }
};