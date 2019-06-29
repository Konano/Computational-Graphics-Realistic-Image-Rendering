/*
* @Author: Konano
* @Date:   2019-06-26 23:16:35
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-26 23:23:37
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <omp.h>
#include <string.h>

int main(int argc, char *argv[])
{
    int w, h, tff;
    FILE *f0 = fopen(argv[1], "r");
    fscanf(f0, "P3\n%d %d\n%d\n", &w, &h, &tff);
    FILE *f1 = fopen(argv[2], "r");
    fscanf(f1, "P3\n%d %d\n%d\n", &w, &h, &tff);
    FILE *f = fopen(argv[3], "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    int a, b, c, tmpa, tmpb, tmpc;
    for (int i = 0; i < w*h; i++) {
        fscanf(f1, "%d%d%d", &a, &b, &c);
        if (a+b+c == 0) fscanf(f0, "%d%d%d", &a, &b, &c); else fscanf(f0, "%d%d%d", &tmpa, &tmpb, &tmpc);
        fprintf(f, "%d %d %d ", a, b, c);
    }
}


