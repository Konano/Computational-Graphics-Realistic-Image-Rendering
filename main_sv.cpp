/*
* @Author: Konano
* @Date:   2019-06-19 15:57:28
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-24 21:55:12
*/

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

    // Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // camera
    // Ray cam(Vec(150, 50, 50), Vec(-1, 0.001, 0.001).norm()); // camera
    Ray cam(Vec(600, 50, 140), Vec(-1, 0.01, -0.1).norm()); // camera
    Vec cx = Vec(0.,0.,-w*.33 / h);
    Vec cy = cx.cross(cam.d).norm()*.33;
    Vec *c = new Vec[w*h]; // image

    // cv::Mat img(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
    // cv::imshow("test", img);

    // system("pause");
    // for (int x = 0; x < h; x++)
    //     for (int y = 0; y < w; y++)
    //         for (int o = 0; o < 3; o++)
    //             img.at<cv::Vec3b>(x, y)[o] = (x < 10 ? 255 : 0);
    // cv::imwrite("output.png", img);
    // return 0;

    Vec r;
    clock_t clock_start = clock(); // Start time
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
// #pragma omp parallel for schedule(dynamic, 1) private(r) num_threads(2) // OpenMP
    for (int y = y0; y < y0+dy; y++) {  // Loop rows
        int sec = (double)(clock() - clock_start) / CLOCKS_PER_SEC / std::max(1,y-y0) * (dy+y0-y);
        fprintf(stderr, "\r%5.2f%%, rows: %d, remain: %02d h %02d m %02d s", 100.*(y-y0)/dy, y, sec/3600, (sec%3600)/60, sec%60);
        // fprintf(stderr, "Thread: %d, Rendering (%d spp) %5.2f%% %d\n", omp_get_thread_num(), samps*4, 100.*(y-y0)/(dy-1), y);
        // fprintf(stderr, "Rendering (%d spp) %5.2f%%\n", samps*4, 100.*y/(h-1));
        // fprintf(stderr, "%d\n", y);
        for (int x = x0; x < x0+dx; x++) {// Loop cols
            // fprintf(stderr, "%d %d\n", y, x);
            for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++) {// 2x2 subpixel rows
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
            // img.at<cv::Vec3b>(h-1-y, x)[2] = gamma(c[(h - y - 1)*w + x].x);
            // img.at<cv::Vec3b>(h-1-y, x)[1] = gamma(c[(h - y - 1)*w + x].y);
            // img.at<cv::Vec3b>(h-1-y, x)[0] = gamma(c[(h - y - 1)*w + x].z);
            // if (omp_get_thread_num() == 0)
            // {
            //     cv::imshow("output", img);
            //     cv::waitKey(1);
            // }
        }
        // cv::imwrite(png_filename, img);
    }
    // printf("\n%d: %d %d\n", NewtonTimes, counter, counter4);
    // for(int i=0; i<NewtonTimes; i++) printf("%03d Times: %d\n", i+1, NewtonCounter[NewtonTimes-i]);
    printf("\n%f sec", (double)(clock() - clock_start) / CLOCKS_PER_SEC);

    FILE *f = fopen(ppm_filename, "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w*h; i++)
        fprintf(f, "%d %d %d ", gamma(c[i].x), gamma(c[i].y), gamma(c[i].z));

    // Save the image.
    // cv::imwrite(png_filename, img);
    // cv::imshow("output", img);
    // cv::waitKey(0);

    // Show the image, wait for user keystroke and quit.
    // cv::imshow("output", img);
    // cv::waitKey(0);
    // cv::destroyAllWindows();
    // system("pause");
}