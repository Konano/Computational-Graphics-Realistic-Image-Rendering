#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include "vec3.hpp"

class Texture {
public:
    Vec color;
    std::string filename;
    unsigned char *buf;
    int w, h, c;
    Texture(const Texture&t): color(t.color), filename(t.filename) {
        if (t.filename != "")
            buf = stbi_load(filename.c_str(), &w, &h, &c, 0);
        else
            buf = NULL;
        // printf("%d %d %d\n", w, h, c);
    }
    Texture(std::string tname, Vec col): color(col), filename(tname) {
        if(tname != "")
            buf = stbi_load(filename.c_str(), &w, &h, &c, 0);
        else
            buf = NULL;
        printf("%d %d %d\n", w, h, c);
    }
    Vec getcol(double a, double b) const {
        if (buf == NULL)
            return color;
        int pw = (int(a * w) % w + w) % w, ph = (int(b * h) % h + h) % h;
        int idx = (h-1-ph) * w * c + pw * c;
        return Vec(buf[idx + 0]+1, buf[idx + 1]+1, buf[idx + 2]+1) / 256.;
    }
};

#endif // __TEXTURE_H__
