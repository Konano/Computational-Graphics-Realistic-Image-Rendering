/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __RENDER_H__
#define __RENDER_H__

#include "scene.hpp"
#include "hitpoint.hpp"

void trace(const Ray &r, Vec fl, int depth, unsigned short* X, HitPoint* hp = NULL) {
    double t;
    int id = 0;
    Vec x, n;
    if (!intersect(r, x, n, id, X)) return;
    const Object *obj = objects[id];
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    Feature ft = feature(obj, x, X);
    Vec f = ft.col;
    double p = f.max();
    if (++depth > 5 || !p) {
        if (erand48(X) < p)
            f = f * (1 / p);
        else
            return;
    }
    switch (ft.rf) {
        case DIFF: {
            if (hp) {
                hp->pos = x;
                hp->fl = fl.mult(f);
                hp->fluxLight += hp->fl.mult(obj->e);
                hp->norm = nl;
                hp->valid = (obj->e.max() == 0);
            } else {
                hitpointsKDTree->update(hitpointsKDTree->root, x, fl, r.d);
                trace(Ray(x, r.d-n*2*n.dot(r.d)), fl.mult(f), depth, X, hp);
            }
            break;
        }
        case SPEC: {
            trace(Ray(x, r.d-n*2*n.dot(r.d)), fl.mult(f), depth, X, hp);
            break;
        }
        case REFR: {
            Ray reflRay(x, r.d - n * 2 * n.dot(r.d));
            bool into = n.dot(nl)>0;
            double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
            if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) {
                trace(reflRay, fl.mult(f), depth, X, hp);
                return;
            }
            Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();

            double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
            double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
            if (depth > 2) {
                if (erand48(X) < P)
                    trace(reflRay, fl.mult(f)*RP, depth, X, hp);
                else
                    trace(Ray(x, tdir), fl.mult(f)*TP, depth, X, hp);
            } else {
                trace(reflRay, fl.mult(f)*Re, depth, X, hp);
                trace(Ray(x, tdir), fl.mult(f)*Tr, depth, X, hp);
            }
            break;
        }
    }
}


// For PT

// Vec radiance(const Ray &r, int depth, unsigned short* X, bool debug = false) {
//     if (debug) puts("radiance.in");
//     if (debug) printf("depth: %d\n", depth);
//     double t;
//     int id = 0;
//     Vec x, n;
//     if (!intersect(r, x, n, id, X))
//         return Vec();
//     const Object *obj = objects[id];
//     if (debug) printf("n x r.d %.12lf\n", n.dot(r.d));
//     Vec nl = n.dot(r.d) < 0 ? n : n * -1;
//     Feature ft = feature(obj, x, X);
//     Vec f = ft.col;
//     if (debug) r.print(), x.print(), n.print(), nl.print(), f.print();
//     double p = f.max();
//     if (debug) puts("pass 28");
//     if (++depth > 5 || !p) {
//         if (erand48(X) < p)
//             f = f * (1 / p);
//         else
//             return obj->e;
//     }
//     if (debug) puts("pass 35");
//     switch (ft.rf) {
//         case DIFF: {
//             double r1 = 2 * M_PI * erand48(X), r2 = erand48(X), r2s = sqrt(r2);
//             Vec w = nl, u = (fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)).cross(w).norm(), v = w.cross(u);
//             Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();
//             if (debug) {
//                 puts("pass DIFF");
//                 Vec ans = obj->e + f.mult(radiance(Ray(x, d), depth, X, debug));
//                 ans.print();
//                 return ans;
//             }
//             return obj->e + f.mult(radiance(Ray(x, d), depth, X, debug));
//         }
//         case SPEC: {
//             if (debug) {
//                 puts("pass SPEC");
//                 Vec ans = obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth, X, debug));
//                 ans.print();
//                 return ans;
//             }
//             return obj->e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth, X, debug));
//         }
//         case REFR: {
//             if (debug) puts("pass REFR");
//             Ray reflRay(x, r.d - n * 2 * n.dot(r.d));
//             if (debug) reflRay.print();
//             bool into = n.dot(nl)>0;
//             double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
//             if (debug) printf("%.6lf %.6lf %.6lf\n", nnt, ddn, 1 - nnt*nnt*(1 - ddn*ddn));
//             if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0)
//                 return obj->e + f.mult(radiance(reflRay, depth, X, debug));
//             Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();

//             double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
//             double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
//             if (debug) printf("%.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", a, b, R0, c, Re, Tr, P, RP, TP);
//             return obj->e + f.mult(depth>2 ? (erand48(X)<P ?
//                 radiance(reflRay, depth, X, debug)*RP : radiance(Ray(x, tdir), depth, X, debug)*TP) :
//                 radiance(reflRay, depth, X, debug)*Re + radiance(Ray(x, tdir), depth, X, debug)*Tr);
//         }
//     }
//     if (debug) puts("radiance.out");
//     return Vec();
// }

#endif // __RENDER_H__