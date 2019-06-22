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
    double dis_to_axis(double l, double r) const {
        if (d.y == 0)
            l = -1e20, r = 1e20;
        else
            l = (l - o.y) / d.y, r = (r - o.y) / d.y;
        if (l > r) std::swap(l, r);

        if (sqr(d.x)+sqr(d.z) == 0)
            return sqr(o.x)+sqr(o.z);
        double p = -(o.x*d.x+o.z*d.z)/(sqr(d.x)+sqr(d.z));
        if (p < l)
            return sqr(pos(l).x)+sqr(pos(l).z);
        if (r < p)
            return sqr(pos(p).x)+sqr(pos(p).z);
        return sqr(o.x)+sqr(o.z)-sqr(o.x*d.x+o.z*d.z)/(sqr(d.x)+sqr(d.z));
    }
    void print() const {
        printf("ray.o: "); o.print();
        printf("ray.d: "); d.print();
    }
};

#endif // __RAY_H__