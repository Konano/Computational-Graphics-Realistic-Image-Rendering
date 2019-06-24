/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __OBJECT_H__
#define __OBJECT_H__

#include "utils.hpp"
#include "vec3.hpp"
#include "ray.hpp"
#include "bezier.hpp"
#include "texture.hpp"

class Object {
public:
    Vec e;     // emission
    Refl_t ty; // material
    Texture texture;

    Object(Vec _e, Vec _c, Refl_t _ty, std::string tname) : e(_e), ty(_ty), texture(tname, _c) {}

    virtual double intersect(const Ray&) const { puts("virtual error in intersect!"); }
    virtual Vec norm(const Vec&) const { puts("virtual error in norm!"); }
    virtual Vec getpos(const Vec&) const { puts("virtual error in getpos!"); }
};

class Sphere : public Object {
public:
    double r;  // radius
    Vec pos;   // position
    Sphere(double _r, Vec _p, Vec _e, Vec _c, Refl_t _ty, std::string tname = "") :
        Object(_e, _c, _ty, tname), r(_r), pos(_p) {}

    virtual double intersect(const Ray& ray) const {
        Vec dis = ray.o - pos;
        double a = ray.d.dot(ray.d);
        double b = (ray.d * 2).dot(dis);
        double c = dis.dot(dis) - r*r;
        double detal =  b * b - 4. * a * c;
        if (detal < 0) return 0;
        double t1 = (- b - sqrt(detal)) / 2 / a;
        double t2 = (- b + sqrt(detal)) / 2 / a;
        if (t1 > eps) return t1;
        else if (t2 > eps) return t2;
        else return 0;
    } // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    virtual Vec norm(const Vec& x) const {
        return (x - pos).norm();
    }
};

class Plane : public Object {
public:
    Vec n, n0;     // store ax+by+cz=1 n=(a,b,c)
    Plane(Vec _n, Vec _e, Vec _c, Refl_t _ty, std::string tname = "") :
        Object(_e, _c, _ty, tname), n(_n), n0(_n.norm()) {}

    virtual double intersect(const Ray &ray) const {
        double t = (1 - ray.o.dot(n)) / ray.d.dot(n);
        if (t < eps) return 0;
        return t;
    }
    virtual Vec norm(const Vec&) const {
        return n0;
    }
};

bool Bezier_debug = true;

class Bezier : public Object {
public:
    BezierCurve2D curve;
    Vec pos; // the buttom center point
    Bezier(Vec _p, BezierCurve2D _c, Vec _e, Vec _col, Refl_t _ty, std::string tname = "") :
        Object(_e, _col, _ty, tname), curve(_c), pos(_p) {}
    Vec getpos(const Vec& o) const {
        double t = curve.dy.solve(o.y - pos.y);
        double u = atan2(o.z - pos.z, o.x - pos.x);
        if (u < 0) u += 2 * PI;
        return Vec(u, t);
    }
    virtual double intersect(const Ray &_ray) const {
        Ray ray = _ray; // copy
        ray.o -= pos;
        double final_t, L, R; bool inside;
        if (ray.in_axis_range(0., curve.height, curve.width, curve.width_mn, L, R, inside) == false) return 0;

        if (fabs(ray.d.y) < 1e-4) {
            // if (Bezier_debug) ray.print();
            // if (Bezier_debug) printf("%.6lf %.6lf\n", fabs(ray.d.y), ray.d.y);
            double hit = ray.o.y;
            for (int i = 0; i <= 3; i++) {
                if (hit < eps || hit > curve.height - eps) return 0;
                // if (Bezier_debug) puts("pass 1");
                double t = curve.dy.solve(hit);
                // if (t < 0 || t > 1) return 0;
                Vec loc = curve.getpos(t);
                double a = sqr(ray.d.x) + sqr(ray.d.z);
                double b = 2. * (ray.o.x * ray.d.x + ray.o.z * ray.d.z);
                double c = sqr(ray.o.x) + sqr(ray.o.z) - sqr(loc.x);
                double detal =  b * b - 4. * a * c;
                if (detal < 0) return 0;
                double t1 = (- b - sqrt(detal)) / 2 / a;
                double t2 = (- b + sqrt(detal)) / 2 / a;
                if (t2 <= eps) return 0;
                // if (Bezier_debug) puts("pass 2");
                if (t1 > eps) final_t = t1; else final_t = t2;
                hit = ray.pos(final_t).y;
            }
            if (Bezier_debug) {
                double tt = curve.dy.solve(hit);
                if (fabs(sqrt(sqr(ray.pos(final_t).x) + sqr(ray.pos(final_t).z)) - curve.dx.calc(tt)) > 1e-6) {
                    // printf("%.12lf\n", fabs(sqrt(sqr(ray.pos(final_t).x) + sqr(ray.pos(final_t).z)) - curve.dx.calc(tt)));
                    // puts("");
                    // ray.print();
                    // printf("final_t: %.12lf hit: %.12lf tt: %.12lf\n", final_t, hit, tt);
                    // ray.pos(final_t).print();
                    // printf("x: %.12lf y: %.12lf\n", curve.dx.calc(tt), curve.dy.calc(tt));
                }
            }
        } else {
            // 常数优化
            // Polynome k = (curve.dy.copy() + (-ray.o.y)) * (1./ray.d.y);
            // Polynome poly = (k.copy() * ray.d.x + ray.o.x).sqr() + (k.copy() * ray.d.z + ray.o.z).sqr() - curve.dx.sqr();
            Polynome k = curve.dy;

            Polynome poly(k.len * 2 - 1);
            for (int i=k.len-1; i>=0; i--) for (int j=k.len-1; j>=0; j--) poly[i+j] -= curve.dx.v[i] * curve.dx.v[j];

            k[0] -= ray.o.y;
            for (int i=k.len-1; i>=0; i--) k[i] /= ray.d.y;
            for (int i=k.len-1; i>=0; i--) k[i] *= ray.d.x;
            k[0] += ray.o.x;
            for (int i=k.len-1; i>=0; i--) for (int j=k.len-1; j>=0; j--) poly[i+j] += k[i] * k[j];

            for (int i=k.len-1; i>=0; i--) k[i] = curve.dy.v[i];
            k[0] -= ray.o.y;
            for (int i=k.len-1; i>=0; i--) k[i] /= ray.d.y;
            for (int i=k.len-1; i>=0; i--) k[i] *= ray.d.z;
            k[0] += ray.o.z;
            for (int i=k.len-1; i>=0; i--) for (int j=k.len-1; j>=0; j--) poly[i+j] += k[i] * k[j];

            // if (Bezier_debug) curve.dy.print();
            // if (Bezier_debug) ray.print();
            // if (Bezier_debug) pos.print();
            // if (Bezier_debug) k.print();
            double ansx;
            double y0 = curve.dy.solve(ray.pos(L).y-eps);
            double y1 = curve.dy.solve(ray.pos(R).y+eps);
            if (inside && (fabs(ray.d.y) < 1e-2)) {
                if (poly.calc(y0-(y1-y0)) < 0 || poly.calc(y1) > 0) {
                    printf("solve_binary %.12lf %.12lf\n", poly.calc(y0-(y1-y0)), poly.calc(y1));
                    return 0;
                } else
                    if (solve_binary(poly, ray, ansx, final_t, y0-(y1-y0), y1) == false) return 0;
            } else {
                if (solve_Newton(poly, ray, ansx, final_t, y0-.1, y1+.1) == false) return 0;
            }

            // for (int i = 1; i <= n-1; i++)
            //     if (ray.dis_to_axis(curve.data[i-1].y0, curve.data[i-1].y1) < curve.data[i-1].width) {

            //     }


            if (Bezier_debug) {
                // puts("pass 161");
                if (fabs(sqrt(sqr(ray.pos(final_t).x) + sqr(ray.pos(final_t).z)) - curve.dx.calc(ansx)) > 1e-2) {
                    // counter3++;
                    printf("%.12lf\n", fabs(sqrt(sqr(ray.pos(final_t).x) + sqr(ray.pos(final_t).z)) - curve.dx.calc(ansx)));
                    // puts("");
                    // ray.print();
                    // printf("final_t: %.12lf y: %.12lf ansx: %.12lf\n", final_t, ray.pos(final_t).y, ansx);
                    // ray.pos(final_t).print();
                    // printf("x: %.12lf y: %.12lf\n", curve.dx.calc(ansx), curve.dy.calc(ansx));
                    // curve.dy.print();
                    // k.print();
                    // poly.print();
                }
            }

        }
        return final_t;
    }

    bool solve_binary(Polynome &poly, Ray &ray, double &ansx, double &t, double l = 0, double r = 1) const {
        if (poly.calc(l) < 0 || poly.calc(r) > 0)
            puts("BUG");
        double mid;
        for (int i = 20; i--;)
            if (poly.calc(mid = (l+r)/2) < 0)
                r = mid;
            else
                l = mid;
        ansx = l;
        t = (curve.dy.calc(ansx)-ray.o.y)/ray.d.y;
        return true;
    }
    bool solve_Newton(Polynome &poly, Ray &ray, double &ansx, double &t, double l = 0, double r = 1) const {
        // if (Bezier_debug) puts("solve_Newton");

        t = 1e20;
        bool fg = false;
        double x = _rand() * (r-l) + l, fx, dx, tmp; //int last = 0;
        // if (poly.constant()) return false;
        for(int i = NewtonTimes, times = 0; i; i--) {
            while (fabs(fx = poly.calc(x)) < eps || fabs(dx = poly.deriv(x)) < eps) {
                // if (Bezier_debug) printf("%.3lf %.3lf %.3lf\n", x, fx, dx);
                if (fabs(fx) < eps) {
                    if ((tmp = (curve.dy.calc(x)-ray.o.y)/ray.d.y) > eps && tmp < t) {
                        // if (fg && t - tmp > 1e-2) counter2++;
                        // if (t - tmp > 1e-2) last = i;
                        ansx = x, t = tmp, (ray.d.y > 0 ? r : l) = x, fg = true;
                    }
                }
                if (++times == 50) return false;
                x = _rand() * (r-l) + l;
            }
            x -= fx / dx;
            if (x < 0 || x > 1) x = _rand() * (r-l) + l;
        }
        // if (fg) counter++; counter4++;
        // if (Bezier_debug && fg) poly.print();
        // if (Bezier_debug && fg) puts("solve_Newton out");
        return fg;
    }
    virtual Vec norm(const Vec& _v) const {
        Vec V = _v - pos;
        double t = curve.dy.solve(V.y);
        double dx = curve.dx.deriv(t);
        double dy = curve.dy.deriv(t);
        double cosv = V.x / sqrt(sqr(V.x) + sqr(V.z));
        if (cosv < -1) cosv = -1; if (cosv > 1) cosv = 1;
        double sinv = sqrt(1 - sqr(cosv)) * (V.z < 0 ? -1 : 1);
        // if (Bezier_debug) {
        //     puts("");
        //     V.print();
        //     printf("%.12f %.12f %.12f\n", dy*cosv, -dx, dy*sinv);
        // }
        return Vec(dy*cosv, -dx, dy*sinv).norm();
    }
};

// class Object {
//     // Texture texture;
//     Object(Texture t): texture(t) {}
//     Object(Refl_t refl, P3 color, P3 emission, ld brdf, std::string tname):
//         texture(tname, brdf, color, emission, refl) {}
//     virtual std::pair<ld, P3> intersect(Ray) {puts("virtual error in intersect!");}
//         // If no intersect, then return (INF, (0,0,0))
//     virtual std::pair<P3, P3> aabb() {puts("virtual error in aabb!");}
//     virtual P3 norm(P3) {puts("virtual error in norm!");}
//         // return norm vec out of obj
//     virtual P3 change_for_bezier(P3) {puts("virtual error in bezier!");}
// };


// struct Sphere
// {
//     double r;  // rad
//     Vec p;     // position
//     Vec e;     // emission
//     Vec c;     // color
//     Refl_t ty; // material

//     Sphere(double _r, Vec _p, Vec _e, Vec _c, Refl_t _ty) : r(_r), p(_p), e(_e), c(_c), ty(_ty) {}

//     double intersect(const Ray &R) const {
//         const double eps = 1e-4;
//         Vec dis = R.o - p;
//         double a = R.d.dot(R.d);
//         double b = (R.d * 2).dot(dis);
//         double c = dis.dot(dis) - r*r;
//         double detal =  b * b - 4. * a * c;
//         if (detal < 0) return 0;
//         double t1 = (- b - sqrt(detal)) / 2 / a;
//         double t2 = (- b + sqrt(detal)) / 2 / a;
//         if (t1 > eps) return t1;
//         else if (t2 > eps) return t2;
//         else return 0;
//     } // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
// };

#endif // __OBJECT_H__