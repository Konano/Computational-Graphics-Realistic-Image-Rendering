/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __FACE_H__
#define __FACE_H__

#include "utils.hpp"
#include "vec3.hpp"

inline bool inTriangle(Vec& a, Vec& b, Vec& c, Vec& p) {
    Vec v0 = c-a;
    Vec v1 = b-a;
    Vec v2 = p-a;

    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    double inverDeno = 1. / (dot00 * dot11 - dot01 * dot01) ;
    double u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
    double v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
    if (u < 0 || u > 1 || v < 0 || v > 1) return false;
    return u + v <= 1;
}

struct Face {
    Vec u, v, w;
    double a, b, c, d;
    Face(Vec _u, Vec _v, Vec _w) : u(_u), v(_v), w(_w) {
        a = (v.y - u.y) * (w.z - u.z) - (w.y - u.y) * (v.z - u.z);
        b = (v.z - u.z) * (w.x - u.x) - (w.z - u.z) * (v.x - u.x);
        c = (v.x - u.x) * (w.y - u.y) - (w.x - u.x) * (v.y - u.y);
        double length = sqrt(sqr(a) + sqr(b) + sqr(c));
        a /= length; b /= length; c /= length;
        d = - a * u.x - b * u.y - c * u.z;
    }
    double intersect(Ray ray, Vec& n) {
        if (a * ray.d.x + b * ray.d.y + c * ray.d.z == 0) return 0;
        double t = (a * ray.o.x + b * ray.o.y + c * ray.o.z + d) / -(a * ray.d.x + b * ray.d.y + c * ray.d.z);
        Vec x = ray.pos(t);
        if (!inTriangle(u, v, w, x)) return 0;
        n = Vec(a, b, c).norm();
        return t;
    }
    Vec min() { return Vmin(Vmin(u, v), w); }
    Vec max() { return Vmax(Vmax(u, v), w); }
};

void loadOBJ(std::string filename, double &ratio, Vec pos, std::vector<Face*>& f) {
    std::vector<Vec> v;
    FILE* file = fopen(filename.c_str(), "r");
    char buf[256];
    double x, y, z;
    int a, b, c;
    Vec tmp = Vec(0,1e100,0);
    while(fscanf(file, "%s", buf) != EOF) {
        switch (buf[0]) {
            case '#':
                fgets(buf, sizeof(buf), file);
                break;
            case 'v':
                fscanf(file, "%lf %lf %lf", &x, &y, &z);
                fgets(buf, sizeof(buf), file);
                v.push_back(Vec((z-x)*sqrt(0.5)*ratio, y*ratio, -(x+z)*sqrt(0.5)*ratio));
                tmp.y = std::min(tmp.y, y*ratio);
                break;
            case 'f':
                fscanf(file, "%d %d %d", &a, &b, &c);
                f.push_back(new Face(v[a-1]+pos-tmp, v[b-1]+pos-tmp, v[c-1]+pos-tmp));
                break;
            default:
                fgets(buf, sizeof(buf), file);
        }
    }
    fclose(file);
}


#endif