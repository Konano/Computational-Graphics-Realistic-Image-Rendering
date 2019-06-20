/*
* @Author: Konano
* @Date:   2019-06-19 15:57:28
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 03:52:19
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define PI 3.14159265358979323846

inline double _rand() { return (double)rand() / (double)RAND_MAX; }
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1.0 / 2.2) * 255 + 0.5); }  // Gamma Correction

struct Vec {
    double x, y, z;
    Vec(double _x = 0, double _y = 0, double _z = 0) : x(_x), y(_y), z(_z) {}
    Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec cross(const Vec &b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    Vec& norm() { return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
};

struct Ray {
    Vec o, d; // origin, direct
    Ray(Vec _o, Vec _d) : o(_o), d(_d) {}
};

enum tyRefl { DIFF, SPEC, GLAS }; //
struct Sphere
{
    double r;  // rad
    Vec p;     // position
    Vec e;     // emission
    Vec c;     // color
    tyRefl ty; // material

    Sphere(double _r, Vec _p, Vec _e, Vec _c, tyRefl _ty) : r(_r), p(_p), e(_e), c(_c), ty(_ty) {}

    double intersect(const Ray &R) const {
        const double eps = 1e-4;
        Vec dis = R.o - p;
        double a = R.d.dot(R.d);
        double b = (R.d * 2).dot(dis);
        double c = dis.dot(dis) - r*r;
        double detal =  b * b - 4. * a * c;
        if (detal < 0) return 0;
        double t1 = (- b - sqrt(detal)) / 2 / a;
        double t2 = (- b + sqrt(detal)) / 2 / a;
        if (t1 > eps) return t1;
        else if (t2 > eps) return t2;
        else return 0;
    } // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
};

Sphere spheres[] = {
    Sphere(1e5,  Vec(1e5+1,40.8,81.6),   Vec(),         Vec(.75,.25,.25), DIFF), //Left
    Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec(),         Vec(.25,.25,.75), DIFF), //Right
    Sphere(1e5,  Vec(50,40.8,1e5),       Vec(),         Vec(.75,.75,.75), DIFF), //Back
    Sphere(1e5,  Vec(50,40.8,-1e5+170),  Vec(),         Vec(),            DIFF), //Front
    Sphere(1e5,  Vec(50,1e5,81.6),       Vec(),         Vec(.75,.75,.75), DIFF), //Bottom
    Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec(),         Vec(.75,.75,.75), DIFF), //Top
    Sphere(16.5, Vec(27,16.5,47),        Vec(),         Vec(1,1,1)*.999,  SPEC), //Mirror
    Sphere(16.5, Vec(73,16.5,78),        Vec(),         Vec(1,1,1)*.999,  GLAS), //Glass
    Sphere(600,  Vec(50,681.6-.27,81.6), Vec(12,12,12), Vec(),            DIFF)  //Light
};

inline bool intersect(const Ray &R, double &t, int &id) {
    int n = (int)sizeof(spheres) / (int)sizeof(Sphere);
    double d, inf = t = 1e20;
    for(int i = n; i--;) {
        if ((d = spheres[i].intersect(R)) && d < t)
            t = d, id = i;
    }
    return t < inf;
}

Vec radiance(const Ray &r, int depth) {
    double t;
    int id = 0;
    if (!intersect(r, t, id))
        return Vec();
    const Sphere &obj = spheres[id];
    Vec x = r.o + r.d*t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1; // nl 与 d 异向
    Vec f = obj.c; // 反射率
    double p = (f.x>f.y && f.x>f.z) ? f.x : (f.y>f.z ? f.y : f.z); // 颜色最大值
    if (++depth > 5 || !p) {
        if (_rand() < p)
            f = f * (1 / p);
        else
            return obj.e;
    }
    switch (obj.ty) {
        case DIFF: {
            double r1 = 2 * M_PI * _rand(), r2 = _rand(), r2s = sqrt(r2);
            Vec w = nl, u = (fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)).cross(w).norm(), v = w.cross(u); // 正交基
            Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();
            return obj.e + f.mult(radiance(Ray(x, d), depth));
        }
        case SPEC: {
            return obj.e + f.mult(radiance(Ray(x, r.d-n*2*n.dot(r.d)), depth));
        }
        case GLAS: {
            Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // 反射光
            bool into = n.dot(nl)>0;  // 光线是否是进入球体
            double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t; // nnt：折射率，ddn：入射角cos
            if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0)    // cos2t<0，无折射
                return obj.e + f.mult(radiance(reflRay, depth));
            Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm(); //折射光

            // 菲涅耳方程
            double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
            double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
            return obj.e + f.mult(depth>2 ? (_rand()<P ?
                radiance(reflRay, depth)*RP : radiance(Ray(x, tdir), depth)*TP) :
                radiance(reflRay, depth)*Re + radiance(Ray(x, tdir), depth)*Tr);
            // RP & TP 不懂
        }
    }
}

int main(int argc, char *argv[]) {
    int w = 1024, h = 768;
    int samps = (argc==2 ? atoi(argv[1])/4 : 1);
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // camera
    Vec cx = Vec(w*.5135 / h)
    Vec cy = cx.cross(cam.d).norm()*.5135;
    Vec *c = new Vec[w*h]; // image

    Vec r;
    clock_t clock_start = clock(); // Start time
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++) {  // Loop rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps*4, 100.*y/(h-1));
        for (unsigned short x = 0; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()) { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {
                        double r1 = 2 * _rand(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * _rand(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0)*(1./samps);
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
                }
    }
    printf("\n%f sec", (double)(clock() - clock_start) / CLOCKS_PER_SEC);

    FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w*h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}