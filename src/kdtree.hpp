/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __KDTREE_H__
#define __KDTREE_H__

#include "face.hpp"

class MeshKDTreeNode {
public:
    Vec min, max;
    std::vector<Face*>* faces;
    MeshKDTreeNode *ls, *rs;
    int l, r;
    bool inside(Face *face) {
        Vec faceMin = face->min();
        Vec faceMax = face->max();
        return (faceMin.x < max.x || faceMin.x == max.x && faceMin.x == faceMax.x)
            && (faceMax.x > min.x || faceMax.x == min.x && faceMin.x == faceMax.x)
            && (faceMin.y < max.y || faceMin.y == max.y && faceMin.y == faceMax.y)
            && (faceMax.y > min.y || faceMax.y == min.y && faceMin.y == faceMax.y)
            && (faceMin.z < max.z || faceMin.z == max.z && faceMin.z == faceMax.z)
            && (faceMax.z > min.z || faceMax.z == min.z && faceMin.z == faceMax.z);
    }
};

class MeshKDTree {
public:
    MeshKDTreeNode* root;
    std::vector<Face*> *faces;
    MeshKDTreeNode* build(int depth, int d, std::vector<Face*>& faces, Vec min, Vec max) {
        MeshKDTreeNode *p = new MeshKDTreeNode;
        p->min = min;
        p->max = max;
        Vec maxL, minR;
        if (d == 0) {
            maxL = Vec((p->min.x + p->max.x) / 2, p->max.y, p->max.z);
            minR = Vec((p->min.x + p->max.x) / 2, p->min.y, p->min.z);
        }
        else if (d == 1) {
            maxL = Vec(p->max.x, (p->min.y + p->max.y) / 2, p->max.z);
            minR = Vec(p->min.x, (p->min.y + p->max.y) / 2, p->min.z);
        }
        else {
            maxL = Vec(p->max.x, p->max.y, (p->min.z + p->max.z) / 2);
            minR = Vec(p->min.x, p->min.y, (p->min.z + p->max.z) / 2);
        }
        p->faces = new std::vector<Face*>;
        for (auto face : faces)
            if (p->inside(face))
                p->faces->push_back(face);

        const int max_faces = 8;
        const int max_depth = 24;

        if (p->faces->size() > max_faces && depth < max_depth) {
            p->ls = build(depth + 1, (d + 1) % 3, *(p->faces), min, maxL);
            p->rs = build(depth + 1, (d + 1) % 3, *(p->faces), minR, max);

            std::vector<Face*> *faceL = p->ls->faces, *faceR = p->rs->faces;
            std::map<Face*, int> cnt;
            for (auto face : *faceL) cnt[face]++;
            for (auto face : *faceR) cnt[face]++;
            p->ls->faces = new std::vector<Face*>;
            p->rs->faces = new std::vector<Face*>;
            p->faces->clear();
            for (auto face : *faceL)
                if (cnt[face] == 1)
                    p->ls->faces->push_back(face);
                else
                    p->faces->push_back(face);
            for (auto face : *faceR)
                if (cnt[face] == 1)
                    p->rs->faces->push_back(face);
        }
        else
            p->ls = p->rs = NULL;
        return p;
    }

    void getFaces(MeshKDTreeNode *p, std::vector<Face*>* faces) {
        p->l = faces->size();
        for (auto face : *(p->faces))
            faces->push_back(face);
        p->r = faces->size();
        if (p->ls) getFaces(p->ls, faces);
        if (p->rs) getFaces(p->rs, faces);
    }

    double getCuboidIntersection(MeshKDTreeNode *p, Ray ray) {
        if (!(ray.o.bigger(p->min) && ray.o.smaller(p->max))) { // outside
            double t = -1e100;
            if (fabs(ray.d.x) > 0)
                t = std::max(t, std::min((p->min.x - ray.o.x) / ray.d.x, (p->max.x - ray.o.x) / ray.d.x));
            if (fabs(ray.d.y) > 0)
                t = std::max(t, std::min((p->min.y - ray.o.y) / ray.d.y, (p->max.y - ray.o.y) / ray.d.y));
            if (fabs(ray.d.z) > 0)
                t = std::max(t, std::min((p->min.z - ray.o.z) / ray.d.z, (p->max.z - ray.o.z) / ray.d.z));
            if (t < eps) return 1e100;
            Vec x = ray.pos(t);
            if (!(x.bigger(p->min) && x.smaller(p->max))) return 1e100;
            return t;
        }
        else return -1e100;
    }

    void getIntersection(MeshKDTreeNode *p, Ray ray, double &tMin, Vec &norm) {
        if (debug)
            counter++;
        int tmpcount = counter;
        if (tmpcount == 7)
            tmpcount = counter;
        Vec n0;
        for (int i = 0; i < p->faces->size(); ++i) {
            double t = (*p->faces)[i]->intersect(ray, n0);
            if (t > 0 && t < tMin) {
                tMin = t;
                norm = n0;
            }
        }

        double tl = p->ls ? getCuboidIntersection(p->ls, ray) : 1e100;
        double tr = p->rs ? getCuboidIntersection(p->rs, ray) : 1e100;
        if (tl < tr) {
            if (tMin <= tl) return;
            if (p->ls) getIntersection(p->ls, ray, tMin, norm);
            if (tMin <= tr) return;
            if (p->rs) getIntersection(p->rs, ray, tMin, norm);
        }
        else {
            if (tMin <= tr) return;
            if (p->rs) getIntersection(p->rs, ray, tMin, norm);
            if (tMin <= tl) return;
            if (p->ls) getIntersection(p->ls, ray, tMin, norm);
        }
    }

    MeshKDTree(std::vector<Face*>& faces) {
        Vec min = Vec(1e100, 1e100, 1e100);
        Vec max = Vec(-1e100, -1e100, -1e100);
        for (auto face : faces) {
            min = Vmin(min, face->min());
            max = Vmax(max, face->max());
        }
        root = build(1, 0, faces, min, max);
        getFaces(root, this->faces = new std::vector<Face*>);
    }
};

#endif // __KDTREE_H__