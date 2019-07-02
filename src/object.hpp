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
#include "kdtree.hpp"

class Object {
public:
    Vec e;     // emission
    Refl_t ty; // material
    Texture texture;

    Object(Vec _e, Vec _c, Refl_t _ty, std::string tname) : e(_e), ty(_ty), texture(tname, _c) {}
    virtual double intersect(const Ray&, Vec &x, Vec &n, unsigned short*) const { puts("virtual error in intersect!"); }
    virtual Vec getpos(const Vec&) const { puts("virtual error in getpos!"); }
};

class Sphere : public Object {
public:
    double r;  // radius
    Vec pos;   // position
    Sphere(double _r, Vec _p, Vec _e, Vec _c, Refl_t _ty, std::string tname = "") :
        Object(_e, _c, _ty, tname), r(_r), pos(_p) {}

    virtual double intersect(const Ray& ray, Vec &x, Vec &n, unsigned short*) const {
        Vec dis = ray.o - pos;
        double a = ray.d.dot(ray.d);
        double b = (ray.d * 2).dot(dis);
        double c = dis.dot(dis) - r*r;
        double detal =  b * b - 4. * a * c;
        if (detal < 0) return 0;
        double t1 = (- b - sqrt(detal)) / 2 / a;
        double t2 = (- b + sqrt(detal)) / 2 / a;
        double t;
        if (t1 > eps) t = t1;
        else if (t2 > eps) t = t2;
        else return 0;
        x = ray.pos(t);
        n = (x - pos).norm();
        return t;
    } // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
};

class Mesh : public Object {
public:
    Vec pos;
    std::vector<Face*> faces;
    MeshKDTree* meshKDTree;
    Mesh(Vec _p, double ratio, Vec _e, Vec _c, Refl_t _ty, std::string obj) :
        Object(_e, _c, _ty, ""), pos(_p) {
            loadOBJ(obj, ratio, pos, faces);
            meshKDTree = new MeshKDTree(faces);
        }

    virtual double intersect(const Ray &ray, Vec &x, Vec &_n, unsigned short*) const {
        if (meshKDTree->getCuboidIntersection(meshKDTree->root, ray) > 1e10)
            return 0;
        double t = 1e100;
        meshKDTree->getIntersection(meshKDTree->root, ray, t, _n);
        x = ray.pos(t);
        return t;
    }
};

class Plane : public Object {
public:
    Vec n, n0;     // store ax+by+cz=1 n=(a,b,c)
    Plane(Vec _n, Vec _e, Vec _c, Refl_t _ty, std::string tname = "") :
        Object(_e, _c, _ty, tname), n(_n), n0(_n.norm()) {}

    virtual double intersect(const Ray &ray, Vec &x, Vec &_n, unsigned short*) const {
        double t = (1 - ray.o.dot(n)) / ray.d.dot(n);
        if (t < eps) return 0;
        x = ray.pos(t);
        _n = n0;
        return t;
    }
};

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
    virtual double intersect(const Ray &_ray, Vec &x, Vec &n, unsigned short *X) const {
        Ray ray = _ray; // copy
        ray.o -= pos;
        double final_t, L, R; bool inside;
        if (ray.in_axis_range(0., curve.height, curve.width, curve.width_mn, L, R, inside) == false) return 0;

        if (fabs(ray.d.y) < 1e-4) {
            double hit = ray.o.y;
            for (int i = 0; i <= 5; i++) {
                if (hit < eps || hit > curve.height - eps) return 0;
                double t = curve.dy.solve(hit);
                Vec loc = curve.getpos(t);
                double a = sqr(ray.d.x) + sqr(ray.d.z);
                double b = 2. * (ray.o.x * ray.d.x + ray.o.z * ray.d.z);
                double c = sqr(ray.o.x) + sqr(ray.o.z) - sqr(loc.x);
                double detal =  b * b - 4. * a * c;
                if (detal < 0) return 0;
                double t1 = (- b - sqrt(detal)) / 2 / a;
                double t2 = (- b + sqrt(detal)) / 2 / a;
                if (t2 <= eps) return 0;
                if (t1 > eps) final_t = t1; else final_t = t2;
                hit = ray.pos(final_t).y;
            }
        } else {
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

            double ansx;
            double y0 = curve.dy.solve(ray.pos(L).y-eps);
            double y1 = curve.dy.solve(ray.pos(R).y+eps);
            if (inside && (fabs(ray.d.y) < 1e-2)) {
                if (poly.calc(y0-(y1-y0)) < 0 || poly.calc(y1) > 0) {
                    // printf("solve_binary %.12lf %.12lf\n", poly.calc(y0-(y1-y0)), poly.calc(y1));
                    return 0;
                } else
                    if (solve_binary(poly, ray, ansx, final_t, y0-(y1-y0), y1) == false) return 0;
            } else {
                if (solve_Newton(poly, ray, ansx, final_t, y0-.1, y1+.1, X) == false) return 0;
            }
        }
        x = _ray.pos(final_t);
        n = norm(x);
        return final_t;
    }

    bool solve_binary(Polynome &poly, Ray &ray, double &ansx, double &t, double l = 0, double r = 1) const {
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
    bool solve_Newton(Polynome &poly, Ray &ray, double &ansx, double &t, double l, double r, unsigned short *X) const {
        l = std::max(l, 0.); r = std::min(r, 1.);
        t = 1e20;
        bool fg = false;
        double x = erand48(X) * (r-l) + l, fx, dx, tmp;
        for(int i = NewtonTimes, times = 0; i>-NewtonTimes; i--, times = 0) {
            while (fabs(fx = poly.calc(x)) < eps || fabs(dx = poly.deriv(x)) < eps) {
                if (fabs(fx) < eps) {
                    if ((tmp = (curve.dy.calc(x)-ray.o.y)/ray.d.y) > eps && tmp < t) {
                        ansx = x, t = tmp, fg = true;
                        (ray.d.y > 0 ? r : l) = x;
                    }
                }
                if (++times == 50)
                    return false;
                x = erand48(X) * (r-l) + l;
            }
            x -= fx / dx;
            if (x < 0 || x > 1) {
                if (i <= 0) break;
                x = erand48(X) * (r-l) + l;
            }
        }
        return fg;
    }
    Vec norm(const Vec& _v) const {
        Vec V = _v - pos;
        double t = curve.dy.solve(V.y);
        double dx = curve.dx.deriv(t);
        double dy = curve.dy.deriv(t);
        double cosv = V.x / sqrt(sqr(V.x) + sqr(V.z));
        if (cosv < -1) cosv = -1; if (cosv > 1) cosv = 1;
        double sinv = sqrt(1 - sqr(cosv)) * (V.z < 0 ? -1 : 1);
        return Vec(dy*cosv, -dx, dy*sinv).norm();
    }
};

#endif // __OBJECT_H__