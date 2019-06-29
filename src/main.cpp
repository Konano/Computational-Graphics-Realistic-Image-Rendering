/*
* @Author: Konano
* @Date:   2019-06-19 15:57:28
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-27 01:09:03
*/

// #define __OPENCV

#include "render.hpp"

// void test() {

//     int h = objects[6]->texture.h;
//     int w = objects[6]->texture.w;

//     cv::Mat img(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
//     for(int y = 0; y < h; y++)
//         for(int x = 0; x < w; x++)
//         {
//             Vec tmp = objects[6]->texture.getcol(1./w*x, 1./h*y);
//             img.at<cv::Vec3b>(y, x)[0] = std::min(int(tmp.x * 255 + 0.5), 255);
//             img.at<cv::Vec3b>(y, x)[1] = std::min(int(tmp.y * 255 + 0.5), 255);
//             img.at<cv::Vec3b>(y, x)[2] = std::min(int(tmp.z * 255 + 0.5), 255);
//         }
//     cv::imshow("output", img);
//     cv::waitKey(0);
// }

int main(int argc, char *argv[]) {

    // printf("%.6lf\n", _rand());

    // test();
    // return 0;

    int samps = (argc>=2 ? atoi(argv[1])/4 : 1);
    char* png_filename = cat(argc>=3 ? argv[2] : "output", ".png");
    char* ppm_filename = cat(argc>=3 ? argv[2] : "output", ".ppm");
    int w = (argc>=5 ? atoi(argv[3]) : 640);
    int h = (argc>=5 ? atoi(argv[4]) : 360);
    int x0 = (argc>=9 ? atoi(argv[5]) : 0);
    int y0 = (argc>=9 ? atoi(argv[6]) : 0);
    int dx = (argc>=9 ? atoi(argv[7]) : w);
    int dy = (argc>=9 ? atoi(argv[8]) : h); y0 = h - (y0+dy);

    // Mirror
    // Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // camera
    // Vec cx = Vec(w*.4 / h);
    // Vec cy = cx.cross(cam.d).norm()*.4;
    // Vec *c = new Vec[w*h]; // image

    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // camera
    Vec cx = Vec(w*.33 / h);
    Vec cy = cx.cross(cam.d).norm()*.33;
    double focusZ = 80; // 焦平面
    double aperture = 1;

    // Ray cam(Vec(600, 50, 140), Vec(-1, 0.01, -0.1).norm()); // camera
    // Vec cx = Vec(0.,0.,-w*.33 / h);
    // Vec cy = cx.cross(cam.d).norm()*.33;

    Vec *c = new Vec[w*h]; // image

#ifdef __OPENCV
    cv::Mat img(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
    // cv::imshow("test", img);
#endif

    // system("pause");
    // for (int x = 0; x < h; x++)
    //     for (int y = 0; y < w; y++)
    //         for (int o = 0; o < 3; o++)
    //             img.at<cv::Vec3b>(x, y)[o] = (x < 10 ? 255 : 0);
    // cv::imwrite("output.png", img);
    // return 0;

    Vec r;
    double clock_start = omp_get_wtime(); // Start time
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
// #pragma omp parallel for schedule(dynamic, 1) private(r) num_threads(2) // OpenMP
    for (int y = y0; y < y0+dy; y++) {  // Loop rows
        int sec = (omp_get_wtime() - clock_start) / std::max(1,y-y0) * (dy+y0-y);
        fprintf(stderr, "\r%5.2f%%, rows: %d, remain: %02d h %02d m %02d s", 100.*(y-y0)/dy, y, sec/3600, (sec%3600)/60, sec%60);
        // fprintf(stderr, "Thread: %d, Rendering (%d spp) %5.2f%% %d\n", omp_get_thread_num(), samps*4, 100.*(y-y0)/(dy-1), y);
        // fprintf(stderr, "Rendering (%d spp) %5.2f%%\n", samps*4, 100.*y/(h-1));
        // fprintf(stderr, "%d\n", y);
        unsigned short X[3] = {0, 0, y*y*y};
        for (int x = x0; x < x0+dx; x++) {// Loop cols
            // printf("%d %d %d %d %d\n", y, x, X[0], X[1], X[2]);
            // if (y == 125 && x == 1244) X[0] = 21280, X[1] = 31644, X[2] = 3421;
            // fprintf(stderr, "%d %d\n", y, x);
            for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++) {// 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()) { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++) {
                        // fprintf(stderr, "\rRendering (%d spp) %5.2f%% %d %d", samps*4, 100.*y/(h-1), y, x);
                         // LINUX
                        // printf("%d %d %d\n", sy, sx, s);
                        double r1 = 2 * erand48(X), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(X), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                            cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        // r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0)*(1./samps);

                        Vec focus = cam.o + d * (focusZ - cam.o.z) / d.z;
                        r1 = 2 * erand48(X); dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        r2 = 2 * erand48(X); dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec cam_o = cam.o + Vec(dx, dy) * aperture;
                        d = (focus - cam_o).norm() * d.len();

                        r = r + radiance(Ray(cam_o + d * 140, d.norm()), 0, X)*(1./samps);
                    }
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
                }
            }
#ifdef __OPENCV
            img.at<cv::Vec3b>(h-1-y, x)[2] = _gamma(c[(h - y - 1)*w + x].x);
            img.at<cv::Vec3b>(h-1-y, x)[1] = _gamma(c[(h - y - 1)*w + x].y);
            img.at<cv::Vec3b>(h-1-y, x)[0] = _gamma(c[(h - y - 1)*w + x].z);
            if (omp_get_thread_num() == 0)
            {
                cv::imshow("output", img);
                cv::waitKey(1);
            }
#endif
        }
        // cv::imwrite(png_filename, img);
    }
    // printf("\n%d: %d %d\n", NewtonTimes, counter, counter4);
    // for(int i=0; i<=NewtonTimes; i++) printf("%03d Times: %d\n", i+1, NewtonCounter[NewtonTimes-i]);
    // for(int i=0; i<=2*NewtonTimes; i++) printf("%03d Times: %d\n", i, NewtonCounter[i]);
    printf("\n%f sec", omp_get_wtime() - clock_start);

    FILE *f = fopen(ppm_filename, "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w*h; i++)
        fprintf(f, "%d %d %d ", _gamma(c[i].x), _gamma(c[i].y), _gamma(c[i].z));

#ifdef __OPENCV
    // Save the image.
    cv::imwrite(png_filename, img);
    // cv::imshow("output", img);
    // cv::waitKey(0);

    // Show the image, wait for user keystroke and quit.
    // cv::imshow("output", img);
    // cv::waitKey(0);
    // cv::destroyAllWindows();
    // system("pause");
#endif
    return 0;
}
