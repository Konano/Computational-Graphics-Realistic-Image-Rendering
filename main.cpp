/*
* @Author: Konano
* @Date:   2019-06-19 15:57:28
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-23 04:40:55
*/

#include "render.hpp"

int main(int argc, char *argv[]) {
    int w = 1024, h = 768;
    int samps = (argc>=2 ? atoi(argv[1])/4 : 1);
    // Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // camera
    Ray cam(Vec(50, 53, 296.6), Vec(0, -0.042612, -1).norm()); // camera
    Vec cx = Vec(w*.5135 / h);
    Vec cy = cx.cross(cam.d).norm()*.5135;
    Vec *c = new Vec[w*h]; // image

    Vec r;
    clock_t clock_start = clock(); // Start time
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
// #pragma omp parallel for schedule(dynamic, 1) private(r) num_threads(5) // OpenMP
    for (int y = 0; y < h; y++) {  // Loop rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps*4, 100.*y/(h-1));
        // fprintf(stderr, "Rendering (%d spp) %5.2f%%\n", samps*4, 100.*y/(h-1));
        for (int x = 0; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()) { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {
                        // fprintf(stderr, "\rRendering (%d spp) %5.2f%% %d %d", samps*4, 100.*y/(h-1), y, x);
                        double r1 = 2 * _rand(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * _rand(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        // r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0)*(1./samps);
                        // if ((y==166 && x==705)) puts("pass main");
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0)*(1./samps);
                        // if ((y==166 && x==705)) puts("out main");
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
                }
    }
    // printf("\n%d: %d %d %d %d\n", NewtonTimes, counter, counter2, counter3, counter4);
    // for(int i=0; i<NewtonTimes; i++) printf("%03d Times: %d\n", i+1, NewtonCounter[NewtonTimes-i]);
    printf("\n%f sec", (double)(clock() - clock_start) / CLOCKS_PER_SEC);

    FILE *f = fopen(argc>=3 ? argv[2] : "image.ppm", "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w*h; i++)
        fprintf(f, "%d %d %d ", gamma(c[i].x), gamma(c[i].y), gamma(c[i].z));
}