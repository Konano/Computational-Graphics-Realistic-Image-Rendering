/*
* @Author: Konano
* @Date:   2019-06-21 13:15:42
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:24:59
*/

#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <string.h>
#include <vector>
#include <map>

#ifdef _WIN32
#define __OPENCV
inline double erand48(unsigned short *X) { return (double)rand() / (double)RAND_MAX; }
#endif

#ifdef __OPENCV
#include "opencv2/opencv.hpp"
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

const double PI = acos(-1);
const double eps = 1e-5;
const double INF = 1e20;

const int checkpoint_interval = 10;

int counter = 0;
// int counter2 = 0;
// int counter3 = 0;
int counter4 = 0;
// int counter5 = 0;
// int counter6 = 0;
int NewtonTimes = 20;
int NewtonCounter[1000];
// bool debug = true;

const double HPR = 100;

const double alpha = 0.7;

double debug = false;

// #define mp(a, b) std::make_pair(a, b)
// typedef std::pair<double,Vec> TN; // t + norm

enum Refl_t { DIFF, SPEC, REFR };

typedef short unsigned su;

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int _gamma(double x) { return std::min(int(pow(clamp(x), 1 / 2.2) * 255 + 0.5), 255); } // Gamma Correction
// inline double gamma_d(double x) { return pow(clamp(x), 1 / 2.2); } // Gamma Correction
inline double sqr(double x) { return x * x; }

#define rep(i, n) for(int i=0; i<n; i++)

inline char* cat(const char* a, const char* b) {
    char *t = new char[strlen(a) + strlen(b) + 1];
    strcpy(t, a);
    strcat(t, b);
    return t;
}

// inline double min3(const double& a, const double& b, const double& c) {
//     return a < b && a < c ? a : (b < c ? b : c);
// }

// inline double max3(const double& a, const double& b, const double& c) {
//     return a > b && a > c ? a : (b > c ? b : c);
// }

#endif // __UTILS_H__
