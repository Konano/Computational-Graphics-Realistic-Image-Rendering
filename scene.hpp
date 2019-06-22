/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __SCENE_H__
#define __SCENE_H__

#include "object.hpp"

#define Bezier_size 9
double Bezier_dx[Bezier_size] = {13.5, 7.5, 10., 12.5, 15., 15., 15., 15., 12.5};
double Bezier_dy[Bezier_size] = {0., 5., 10., 20., 25., 30., 35., 40., 40.3};

int objects_num = 8;
Object* objects[] = {
// Object objects[] = {
    new Plane(Vec(1,     0,       0),      Vec(), Vec(.75,.25,.25), DIFF),
    new Plane(Vec(1./99, 0,       0),      Vec(), Vec(.25,.25,.75), DIFF),
    new Plane(Vec(0,     0,       1),      Vec(), Vec(.75,.75,.75), DIFF),
    new Plane(Vec(0,     0,       1./171), Vec(), Vec(),            DIFF),
    new Plane(Vec(0,     1,       0),      Vec(), Vec(.75,.75,.75), DIFF),
    new Plane(Vec(0,     1./82.6, 0),      Vec(), Vec(.75,.75,.75), DIFF),
    // new Sphere(16.5, Vec(27, 17.5, 48),        Vec(),         Vec(1,1,1)*.999,SPEC), //Mirror
    // new Sphere(16.5, Vec(73, 17.5, 79),        Vec(),         Vec(1,1,1)*.999,SPEC),
    new Bezier(Vec(73, 0, 79), BezierCurve2D(Bezier_dx, Bezier_dy, Bezier_size), Vec(), Vec(1,1,1)*.999, SPEC),
    // new Sphere(600,  Vec(50, 682.6-.27, 82.6), Vec(12,12,12), Vec(),          DIFF)  //Light
    new Sphere(1.5,  Vec(50, 82.6-16.5, 82.6), Vec(4,4,4)*100,Vec(),          DIFF),//Lite

    // new Sphere(1e5,  Vec(1e5+1,40.8,81.6),   Vec(),         Vec(.75,.25,.25), DIFF), //Left
    // new Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec(),         Vec(.25,.25,.75), DIFF), //Right
    // new Sphere(1e5,  Vec(50,40.8,1e5),       Vec(),         Vec(.75,.75,.75), DIFF), //Back
    // new Sphere(1e5,  Vec(50,40.8,-1e5+170),  Vec(),         Vec(),            DIFF), //Front
    // new Sphere(1e5,  Vec(50,1e5,81.6),       Vec(),         Vec(.75,.75,.75), DIFF), //Bottom
    // new Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec(),         Vec(.75,.75,.75), DIFF), //Top
    // new Sphere(16.5, Vec(27, 16.5, 47),        Vec(),         Vec(1,1,1)*.999,SPEC), //Mirror
    // new Sphere(16.5, Vec(73, 16.5, 78),        Vec(),         Vec(1,1,1)*.999,REFR), //REFRs
    // new Sphere(600,  Vec(50, 681.6-.27, 81.6), Vec(12,12,12), Vec(),          DIFF)  //Light
};

inline bool intersect(const Ray &R, double &t, int &id, bool debug = false) {
    double d, inf = t = 1e20;
    for(int i = objects_num; i--;) {
        if (debug) puts("intersect.in"), Bezier_debug = true;
        d = objects[i]->intersect(R);
        if (debug) puts("intersect.out"), Bezier_debug = false;
        if (d && d < t)
        // if ((d = objects[i]->intersect(R)) && d < t)
            t = d, id = i;
    }
    // if (debug) puts("Exit");
    return t < inf;
}

#endif // __SCENE_H__