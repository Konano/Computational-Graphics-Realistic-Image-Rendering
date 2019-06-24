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
#include <string.h>
#include <opencv2/opencv.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

const double PI = acos(-1);
const double eps = 1e-5;
const double INF = 1e20;

int counter = 0;
// int counter2 = 0;
// int counter3 = 0;
int counter4 = 0;
// int counter5 = 0;
// int counter6 = 0;
int NewtonTimes = 25;
// int NewtonCounter[1000];
// bool debug = true;

// #define mp(a, b) std::make_pair(a, b)
// typedef std::pair<double,Vec> TN; // t + norm

enum Refl_t { DIFF, SPEC, REFR };

inline double _rand() { return (double)rand() / (double)RAND_MAX; }
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int gamma(double x) { return std::min(int(pow(clamp(x), 1 / 2.2) * 255 + 0.5), 255); } // Gamma Correction
// inline double gamma_d(double x) { return pow(clamp(x), 1 / 2.2); } // Gamma Correction
inline double sqr(double x) { return x * x; }

inline char* cat(const char* a, const char* b) {
    char *t = new char[strlen(a) + strlen(b) + 1];
    strcpy(t, a);
    strcat(t, b);
    return t;
}

#endif // __UTILS_H__
