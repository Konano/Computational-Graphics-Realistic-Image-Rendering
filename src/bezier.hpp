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
    double height, width, width_mn;
    int n;
    BezierCurve2D(double* px, double* py, int _n): dx(Polynome(_n)), dy(Polynome(_n)), n(_n) {
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
        width_mn = 1e20;
        for (double t = 0; t <= 1; t += 0.00001) {
            Vec pos = getpos(t);
            if (width < pos.x) width = pos.x;
            if (width_mn > pos.x) width_mn = pos.x;
        }
    }
    Vec getpos(double t) const { return Vec(dx.calc(t), dy.calc(t)); }
};

#endif // __BEZIER_H__