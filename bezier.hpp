/*
* @Author: Konano
* @Date:   2019-06-21 13:31:34
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:31:34
*/

#ifndef __BEZIER_H__
#define __BEZIER_H__

#include "utils.hpp"
#include "polynome.hpp"

class BezierCurve2D {
public:
    Polynome dx, dy;
    double height, width;
    int n;
    // int num; double r;
    // struct D {
    //     double t0, t1, y0, y1, width;
    // } data[10];
    BezierCurve2D(double* px, double* py, int _n): dx(Polynome(_n)), dy(Polynome(_n)), n(_n) {
    // BezierCurve2D(double* px, double* py, int _n, int _num, double _r): num(num_), n(n_), r(r_) {
        n--;
        for(int i = 0, C = 1; i <= n; i++) {
            C = C * (i ? n-i+1 : 1) / (i ? i : 1);
            for(int j = 0, C0 = 1; j <= n-i; j++) {
                C0 = C0 * (j ? n-i-j+1 : 1) / (j ? j : 1);
                dx[j+i] += px[i] * C * C0 * (j&1 ? -1 : 1);
                dy[j+i] += py[i] * C * C0 * (j&1 ? -1 : 1);
            }
        }
        height = dy.calc(1);

        width = 0;
        for (double t = 0; t <= 1; t += 0.00001) {
            Vec pos = getpos(t);
            if (width < pos.x) width = pos.x;
        }

        // width = 0;
        // double interval = 1. / n, c = interval;
        // for (int i = 1; i <= n-1; c += interval, i++) {
        //     D &dt = data[i-1];
        //     dt.width = 0;
        //     dt.t0 = std::max(0., c - interval);
        //     dt.t1 = std::min(1., c + interval);
        //     dt.y0 = getpos(dt.t0).y;
        //     dt.y1 = getpos(dt.t1).y;
        //     for (double t = dt.t0; t <= dt.t1; t += 0.00001) {
        //         Vec pos = getpos(t);
        //         if (dt.width < pos.x) dt.width = pos.x;
        //     }
        //     dt.width += eps;
        //     if (width < dt.width) width = dt.width;
        // }
    }
    Vec getpos(double t) const { return Vec(dx.calc(t), dy.calc(t)); }
};

#endif // __BEZIER_H__