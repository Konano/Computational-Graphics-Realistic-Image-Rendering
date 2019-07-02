/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __SCENE_H__
#define __SCENE_H__

#include "object.hpp"

// #define Bezier_size 9
// double Bezier_dx[Bezier_size] = {13.5, 7.5, 10., 12.5, 15., 15., 15., 15., 12.5};
// double Bezier_dy[Bezier_size] = {0., 5., 10., 20., 25., 30., 35., 40., 40.3};

struct Feature {
public:
    Vec col;
    Refl_t rf;
    Feature(Vec _col, Refl_t _rf) : col(_col), rf(_rf) {}
};

#define Bezier_size 4
double Bezier_dx[Bezier_size] = {15, 42, -6, 12};
double Bezier_dy[Bezier_size] = {0, 24, 46.2, 72.6};
// double Bezier_dx[Bezier_size] = {5.5, 5.5, 5.5, 5.5};
// double Bezier_dy[Bezier_size] = {0, 20, 38.5, 55.5};

// int objects_num = 8;
int objects_num = 14;
// int objects_num = 10;
Object* objects[] = {
// Object objects[] = {


    new Plane(Vec(1,     0,       0),      Vec(), Vec(.75,.25,.25), DIFF, "wall_front.png"),
    new Plane(Vec(1./478, 0,       0),      Vec(), Vec(), DIFF),
    new Plane(Vec(0,     0,       1),      Vec(), Vec(.25,.25,.75), DIFF, "wall_right.jpg"),
    new Plane(Vec(0,     0,       1./293), Vec(), Vec(.8,.8,.8), DIFF),
    new Plane(Vec(0,     1,       0),      Vec(), Vec(.25,.75,.25), DIFF, "floor.png"),
    new Plane(Vec(0,     1./200, 0),      Vec(), Vec(.8,.8,.8), DIFF),
    new Sphere(14, Vec(196, 15+.0001, 23),        Vec(),         Vec(.5,.5,1),REFR),
    new Sphere(8, Vec(266, 9+.0001, 70),        Vec(),         Vec(.5,1,.5),REFR),
    new Sphere(8, Vec(34, 9+.0001, 166+30),        Vec(),         Vec(1,.5,.5),REFR),
    // new Sphere(8, Vec(145, 9+.0001, 84),        Vec(),         Vec(.8,.4,.8),REFR),
    new Sphere(14, Vec(111, 15+.0001, 132+25),        Vec(),         Vec(1,1,1)*.999,DIFF,"cueball.png"),
    new Sphere(25, Vec(86, 26+.0001, 97+30),        Vec(),         Vec(1,1,1)*.999,DIFF,"blueball.jpg"),
    // new Sphere(16.5, Vec(73, 17.5, 79),        Vec(),         Vec(1,1,1)*.999,SPEC),
    new Bezier(Vec(35, 1, 30), BezierCurve2D(Bezier_dx, Bezier_dy, Bezier_size), Vec(), Vec(.75, .25, .25), DIFF, "vase.png"),
    new Sphere(1200,  Vec(160, 1400-.27, 120), Vec(1,1,1), Vec(1,1,1)*.8,          DIFF),  //Light
    // new Sphere(6,  Vec(100, 150, 100), Vec(4,4,4)*100, Vec(),          DIFF),//Lite
    new Mesh(Vec(120, 1+.0001, 84), 0.1, Vec(), Vec(1,1,1)*.9, DIFF, "angel_lucy.obj")


    // new Sphere(1e5,  Vec(1e5+1,40.8,81.6),   Vec(),         Vec(1,.9,1)*.7, SPEC), //Left
    // new Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec(),         Vec(.9,1,.9)*.7, SPEC), //Right
    // new Sphere(1e5,  Vec(50,40.8,1e5),       Vec(),         Vec(.9,1,1)*.7, SPEC), //Back
    // new Sphere(1e5,  Vec(50,40.8,-1e5+170),  Vec(),         Vec(1,1,.9)*.7, SPEC), //Front
    // new Sphere(1e5,  Vec(50,1e5,81.6),       Vec(),         Vec(1,.9,.9)*.7, SPEC), //Bottom
    // new Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec(),         Vec(.9,.9,1)*.7, SPEC), //Top
    // new Sphere(12, Vec(27, 65, 47),        Vec(.03,0,0),         Vec(1,1,1)*.999  ,REFR),
    // new Sphere(9, Vec(12, 30, 37),        Vec(0,.05,0),          Vec(1,1,1)*.999  ,REFR),
    // new Sphere(6, Vec(73, 36, 78),        Vec(0,0,.07),          Vec(1,1,1)*.999   ,REFR), //REFRs
    // new Sphere(10, Vec(50, 12, 80),        Vec(),         Vec(1,1,1)*.999   ,REFR), //REFRs

    // new Plane(Vec(1,     0,       0),      Vec(), Vec(.75,.25,.25), DIFF, "wall_front.png"),
    // new Plane(Vec(1./478, 0,       0),      Vec(), Vec(), DIFF),
    // new Plane(Vec(0,     0,       1),      Vec(), Vec(.25,.25,.75), DIFF, "wall_right.jpg"),
    // new Plane(Vec(0,     0,       1./293), Vec(), Vec(.8,.8,.8), DIFF),
    // new Plane(Vec(0,     1,       0),      Vec(), Vec(.25,.75,.25), DIFF, "floor.png"),
    // new Plane(Vec(0,     1./200, 0),      Vec(), Vec(.8,.8,.8), DIFF),
    // new Sphere(1200,  Vec(145, 1400-.27, 99), Vec(30,30,30), Vec(), DIFF),  //Light
    // new Mesh(Vec(145, 1+.0001, 99), 0.1, Vec(), Vec(1,1,1)*.9, DIFF, "angel_lucy.obj")
};

const double LightR = 25.454412191209601766461707941066;
const double LightX = 160;
const double LightY = 200;
const double LightZ = 120;

inline bool intersect(const Ray &R, Vec &x, Vec &n, int &id, unsigned short* X) {
    double d, inf = 1e20, t = 1e20; Vec x0, n0;
    for(int i = objects_num; i--;) {
        d = objects[i]->intersect(R, x0, n0, X);
        if (d && d < t)
            t = d, id = i, x = x0, n = n0;
    }
    return t < inf;
}

Feature feature(const Object* obj, Vec o, unsigned short* X) {
    const Texture &texture = obj->texture;
    if (texture.filename == "vase.png") {
        Vec tmp = obj->getpos(o);
        if (erand48(X) < 0.05)
            return Feature(Vec(1,1,1)*.99,
                           SPEC);
        return Feature(texture.getcol(tmp.x/2/PI+.5, tmp.y),
                       obj->ty);
    }
    if (texture.filename == "wall_right.jpg") {
        return Feature(texture.getcol(-o.x*12/3143, -o.y*12/1834),
                       obj->ty);
    }
    if (texture.filename == "floor.png") {
        Vec col = texture.getcol(o.x / 104, o.z / 60);
        if (int(col.x*256-1) == 233 && int(col.y*256-1) == 233 && int(col.z*256-1) == 233)
            return Feature(Vec(1,1,1)*.999,
                           SPEC);
        return Feature(col,
                       obj->ty);
    }
    if (texture.filename == "wall_front.png") {
        return Feature(texture.getcol(o.z*6/1920, o.y*6/1200),
                       obj->ty);
    }
    if (texture.filename == "cueball.png") {
        o = (o - Vec(111, 15, 132+25)) / 14;
        o.norm();
        Vec x = Vec(11, 45, 14).norm(); // e chou
        Vec y = x.cross(Vec(19, 19, 81).norm()).norm();
        Vec z = y.cross(x).norm();
        if (erand48(X) < 0.02)
            return Feature(Vec(1,1,1)*.999,
                           SPEC);
        return Feature(texture.getcol(o.dot(x)*.5+.5, o.dot(z)*.5+.5),
                       obj->ty);
    }
    if (texture.filename == "ball.jpg") {
        o = (o - Vec(86, 26, 97)) / 25;
        o.norm();
        Vec x = Vec(114, 51, 41).norm(); // e chou
        Vec y = x.cross(Vec(91, 91, 80).norm()).norm();
        Vec z = y.cross(x).norm();
        return Feature(texture.getcol(acos(o.dot(z))/PI+0.5, atan(o.dot(y)/o.dot(x))/2/PI),
                       obj->ty);
    }
    if (texture.filename == "blueball.jpg") {
        o = (o - Vec(86, 26, 97+30)) / 25;
        o.norm();
        Vec x = Vec(-0.2, 1, 0.3).norm();
        Vec y = x.cross(Vec(91, 91, 80).norm()).norm();
        Vec z = y.cross(x).norm();
        double fx = acos(o.dot(x))/PI;
        double fy = atan(o.dot(y)/o.dot(z))/2/PI+0.5;
        int dx = floor(fx*16);
        int dy = floor(fy*32);
        int o = floor(dx+1 - fx*16 + fy*32 - dy)+1;
        dy %= 16;
        unsigned short XX[3] = {(su)(dy*dy*o), (su)(dx*dx*o), (su)(dx*dy*o)};
        int colid = floor(erand48(XX)*5); Vec col;
        switch (colid) {
            case 0: col = Vec(45, 125, 185); break;
            case 1: col = Vec(95, 155, 200); break;
            case 2: col = Vec(145, 185, 215); break;
            case 3: col = Vec(195, 215, 230); break;
            case 4: col = Vec(245, 245, 245); break;
        }
        return Feature(col * (1./255),
                       obj->ty);
    }

    return Feature(texture.color,
                   obj->ty);
}

inline Ray generateRay(unsigned short* X) {

    double dx = 2 * erand48(X) - 1, dz = 2 * erand48(X) - 1;
    while (sqr(dx) + sqr(dz) > 1) dx = 2 * erand48(X) - 1, dz = 2 * erand48(X) - 1;

    double rd = erand48(X) * PI, rd2 = erand48(X) * PI * 2;
    Vec o = Vec(LightX + dx*LightR, LightY, LightZ + dz*LightR);
    Vec d = Vec(cos(rd)*cos(rd2), -sin(rd), cos(rd)*sin(rd2));

    return Ray(o + d, d);
}

#endif // __SCENE_H__