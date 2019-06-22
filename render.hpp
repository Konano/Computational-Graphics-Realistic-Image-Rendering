/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __RENDER_H__
#define __RENDER_H__

#include "scene.hpp"

Vec radiance(const Ray &r, int depth, bool debug = false) {
    if (debug) puts("radiance.in");
    if (debug) printf("depth: %d\n", depth);
    double t;
    int id = 0;
    if (!intersect(r, t, id, (depth==8)&&debug))
        return Vec();
    const Object *obj = objects[id];
    Vec x = r.o + r.d*t;
    // Vec n = (x - obj.p).norm();
    Vec n = obj->norm(x);
    Vec nl = n.dot(r.d) < 0 ? n : n * -1; // nl 与 d 异向
    Vec f = obj->c; // 反射率
    double p = f.max(); // 颜色最大值
    if (debug) puts("pass 28");
    if (++depth > 5 || !p) {
        if (_rand() < p)
            f = f * (1 / p);
        else
            return obj->e;
    }
    if (debug) puts("pass 35");
    switch (obj->ty) {
        case DIFF: {
            if (debug) puts("pass DIFF");
            double r1 = 2 * M_PI * _rand(), r2 = _rand(), r2s = sqrt(r2);
            Vec w = nl, u = (fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)).cross(w).norm(), v = w.cross(u); // 正交基
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();
            return obj->e + f.mult(radiance(Ray(x, d), depth, debug));
        }
        case SPEC: {
            if (debug) puts("pass SPEC");
            // return obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth));
            return obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth, debug));
        }
        case REFR: {
            puts("pass REFR");
            Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // 反射光
            bool into = n.dot(nl)>0;  // 光线是否是进入球体
            double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t; // nnt：折射率，ddn：入射角cos
            if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0)    // cos2t<0，无折射
                return obj->e + f.mult(radiance(reflRay, depth));
            Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm(); //折射光

            // 菲涅耳方程
            double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
            double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
            return obj->e + f.mult(depth>2 ? (_rand()<P ?
                radiance(reflRay, depth)*RP : radiance(Ray(x, tdir), depth)*TP) :
                radiance(reflRay, depth)*Re + radiance(Ray(x, tdir), depth)*Tr); // 俄罗斯轮盘赌
        }
    }
    if (debug) puts("radiance.out");
    return Vec();
}

#endif // __RENDER_H__