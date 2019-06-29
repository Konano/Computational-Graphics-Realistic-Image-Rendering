/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __RENDER_H__
#define __RENDER_H__

#include "scene.hpp"

Vec radiance(const Ray &r, int depth, unsigned short* X, bool debug = false) {
    if (debug) puts("radiance.in");
    if (debug) printf("depth: %d\n", depth);
    double t;
    int id = 0;
    if (!intersect(r, t, id, X))
        return Vec();
    if (debug) printf("t: %.6lf, id: %d\n", t, id);
    const Object *obj = objects[id];
    Vec x = r.o + r.d*t;
    // Vec n = (x - obj.p).norm();
    Vec n = obj->norm(x);
    if (debug) printf("n x r.d %.12lf\n", n.dot(r.d));
    Vec nl = n.dot(r.d) < 0 ? n : n * -1; // nl 与 d 异向
    Feature ft = feature(obj, x, X);
    Vec f = ft.col; // 反射率
    // if (id == 6) f.print();
    if (debug) r.print(), x.print(), n.print(), nl.print(), f.print();
    double p = f.max(); // 颜色最大值
    if (debug) puts("pass 28");
    if (++depth > 5 || !p) {
        if (erand48(X) < p)
            f = f * (1 / p);
        else
            // return obj->e;
            {if (debug) obj->e.print();return obj->e;}
    }
    if (debug) puts("pass 35");
    switch (ft.rf) {
        case DIFF: {
            double r1 = 2 * M_PI * erand48(X), r2 = erand48(X), r2s = sqrt(r2);
            Vec w = nl, u = (fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)).cross(w).norm(), v = w.cross(u); // 正交基
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();
            if (debug) {
                puts("pass DIFF");
                Vec ans = obj->e + f.mult(radiance(Ray(x, d), depth, X, debug));
                ans.print();
                return ans;
            }
            return obj->e + f.mult(radiance(Ray(x, d), depth, X, debug));
        }
        case SPEC: {
            if (debug) {
                puts("pass SPEC");
                Vec ans = obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth, X, debug));
                ans.print();
                return ans;
            }
            // return obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth));
            return obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth, X, debug));
        }
        case REFR: {
            if (debug) puts("pass REFR");
            Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // 反射光
            if (debug) reflRay.print();
            bool into = n.dot(nl)>0;  // 光线是否是进入球体
            double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t; // nnt：折射率，ddn：入射角cos
            if (debug) printf("%.6lf %.6lf %.6lf\n", nnt, ddn, 1 - nnt*nnt*(1 - ddn*ddn));
            if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0)    // cos2t<0，无折射
                return obj->e + f.mult(radiance(reflRay, depth, X, debug));
            Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm(); //折射光

            // 菲涅耳方程
            double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
            double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
            if (debug) printf("%.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", a, b, R0, c, Re, Tr, P, RP, TP);
            return obj->e + f.mult(depth>2 ? (erand48(X)<P ?
                radiance(reflRay, depth, X, debug)*RP : radiance(Ray(x, tdir), depth, X, debug)*TP) :
                radiance(reflRay, depth, X, debug)*Re + radiance(Ray(x, tdir), depth, X, debug)*Tr); // 俄罗斯轮盘赌
        }
    }
    if (debug) puts("radiance.out");
    return Vec();
}

#endif // __RENDER_H__