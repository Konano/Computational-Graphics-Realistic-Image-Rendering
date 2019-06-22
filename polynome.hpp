/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __POLYNOME_H__
#define __POLYNOME_H__

#include "utils.hpp"

class Polynome
{
public:
    int len;
    double *v;
    Polynome() { len = 0; v = NULL; }
    Polynome(int _l) : len(_l) {
        v = new double[len];
        for(int i = 0; i < len; i++) v[i] = 0;
        // time_count += len;
    }
    Polynome(const Polynome& p)
    {
        len = p.len;
        v = new double[len];
        for(int i = 0; i < len; i++) v[i] = p.v[i];
        // time_count += len;
    }
    ~Polynome() { if (v != NULL) delete [] v; }

    double& operator[](int i) { return v[i]; }

    Polynome copy() const { return Polynome(*this); }
    Polynome operator+(double a) {
        v[0] += a;
        return *this;
    }
    Polynome operator*(double a) {
        for(int i = 0; i < len; i++) v[i] *= a;
        return *this;
    }
    Polynome operator+(const Polynome& b) {
        for(int i = 0; i < len; i++) v[i] += b.v[i];
        return *this;
    }
    Polynome operator-(const Polynome& b) {
        for(int i = 0; i < len; i++) v[i] -= b.v[i];
        return *this;
    }
    Polynome sqr() const {
        Polynome tmp = Polynome(len * 2 - 1);
        for (int i = 0; i < len; i++)
            for (int j = 0; j < len; j++)
                tmp[i+j] += v[i] * v[j];
        return tmp;
    }

    double calc(double x) const {
        double sum = 0, pow = 1;
        for(int i = 0; i < len; i++) {
            sum += pow * v[i];
            pow *= x;
        }
        return sum;
    }
    double deriv(double x) const {
        double sum = 0, pow = 1;
        for(int i = 1; i < len; i++) {
            sum += pow * v[i] * i;
            pow *= x;
        }
        return sum;
    }
    double solve(double v, double l = 0, double r = 1) const {
        // Sepcial for y
        double mid;
        for (int i=20; i; i--) {
            mid = (l+r) / 2;
            if (calc(mid) < v) l = mid; else r = mid;
        }
        return l;
    }
    void print() const {
        for(int i = 0; i < len; i++) printf("%.12lf%c", v[i], i==len-1?'\n':' ');
    }
};

#endif