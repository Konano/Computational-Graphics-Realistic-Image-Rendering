/*
* @Author: Konano
* @Date:   2019-06-21 13:36:08
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:36:08
*/

#ifndef __RAY_H__
#define __RAY_H__


#include "vec3.hpp"

struct Ray {
    Vec o, d; // origin, direct
    Ray(Vec _o, Vec _d) : o(_o), d(_d) {}
    Vec pos(double t) const { return o + d*t; }
    bool in_axis_range(double l, double r, double radmx, double radmn, double &L, double &R, bool &inside) const {
        if (d.y == 0)
            L = 0, R = 1e20;
        else
            L = (l - o.y) / d.y, R = (r - o.y) / d.y;
        if (L > R) std::swap(L, R);
        if (L < 0) L = 0;

        if (sqr(d.x)+sqr(d.z) == 0)
            return sqr(o.x)+sqr(o.z) < sqr(radmx);

        double dis;

        double p = -(o.x*d.x+o.z*d.z)/(sqr(d.x)+sqr(d.z));
        if (p < L)
            dis = sqr(pos(L).x)+sqr(pos(L).z);
        else if (R < p)
            dis = sqr(pos(R).x)+sqr(pos(R).z);
        else
            dis = sqr(o.x)+sqr(o.z)-sqr(o.x*d.x+o.z*d.z)/(sqr(d.x)+sqr(d.z));
        if (dis > sqr(radmx)) return false;

        double L0 = p - sqrt(p*p - (sqr(o.x)+sqr(o.z)-sqr(radmx))/(sqr(d.x)+sqr(d.z)));
        double R0 = p + sqrt(p*p - (sqr(o.x)+sqr(o.z)-sqr(radmx))/(sqr(d.x)+sqr(d.z)));
        if (L < L0 - eps && dis < sqr(radmn)) {
            L = std::max(L, L0);
            R = std::min(R, p);
            inside = true;
        } else {
            L = std::max(L, L0);
            R = std::min(R, R0);
            inside = false;
        }

        return true;
    }
    void print() const {
        printf("ray.o: "); o.print();
        printf("ray.d: "); d.print();
    }
};

#endif // __RAY_H__