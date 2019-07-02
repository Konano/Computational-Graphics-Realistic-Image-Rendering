/*
* @Author: Konano
* @Date:   2019-06-21 13:25:19
* @Last Modified by:   Konano
* @Last Modified time: 2019-06-21 13:27:04
*/

#ifndef __HITPOINT_H__
#define __HITPOINT_H__

struct HitPoint {
    Vec pos, fl, flux, fluxLight, norm;
    double r2;
    int n;
    bool valid;
    HitPoint() {
        pos = fl = flux = fluxLight = norm = Vec();
        valid = false; n = 0;
        r2 = HPR;
    }
};

bool cmpHitPointX(HitPoint *a, HitPoint *b) {
    return a->pos.x < b->pos.x;
}

bool cmpHitPointY(HitPoint *a, HitPoint *b) {
    return a->pos.y < b->pos.y;
}

bool cmpHitPointZ(HitPoint *a, HitPoint *b) {
    return a->pos.z < b->pos.z;
}

struct HitPointKDTreeNode {
    HitPoint *hitpoint;
    Vec min, max;
    double maxr2;
    HitPointKDTreeNode *ls, *rs;
};

class HitPointKDTree {
private:
    int n;
    HitPoint** hitpoints;
    void build(HitPointKDTreeNode *p, int l, int r, int d) {
        p->min = Vec(1e100, 1e100, 1e100);
        p->max = Vec(-1e100, -1e100, -1e100);
        p->maxr2 = 0;
        for (int i = l; i <= r; ++i) {
            p->min = Vmin(p->min, hitpoints[i]->pos);
            p->max = Vmax(p->max, hitpoints[i]->pos);
            p->maxr2 = std::max(p->maxr2, hitpoints[i]->r2);
        }
        int m = (l+r)>>1;
        if (d == 0)
            std::nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointX);
        else if (d == 1)
            std::nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointY);
        else
            std::nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointZ);
        p->hitpoint = hitpoints[m];
        if (l <= m - 1)
            #pragma omp task
            build(p->ls = new HitPointKDTreeNode, l, m-1, (d+1)%3);
        else
            p->ls = NULL;
        if (m + 1 <= r)
            #pragma omp task
            build(p->rs = new HitPointKDTreeNode, m+1, r, (d+1)%3);
        else
            p->rs = NULL;
    }
    void del(HitPointKDTreeNode *p) {
        if (p->ls) del(p->ls);
        if (p->rs) del(p->rs);
        delete p;
    }
public:
    HitPointKDTreeNode *root;
    HitPointKDTree(std::vector<HitPoint*>& hitpoints) {
        n = hitpoints.size();
        this->hitpoints = new HitPoint*[n];
        rep(i, n) this->hitpoints[i] = hitpoints[i];
        #pragma omp parallel
        {
            #pragma omp single nowait
            build(root = new HitPointKDTreeNode, 0, n-1, 0);
            #pragma omp barrier
        }
    }
    ~HitPointKDTree() {
        if (!root) return;
        del(root);
        delete[] hitpoints;
    }
    void update(HitPointKDTreeNode *p, Vec pos, Vec fl, Vec d) {
        if (!p) return;
        double mind = 0;
        if (pos.x > p->max.x) mind += sqr(pos.x - p->max.x);
        if (pos.x < p->min.x) mind += sqr(p->min.x - pos.x);
        if (pos.y > p->max.y) mind += sqr(pos.y - p->max.y);
        if (pos.y < p->min.y) mind += sqr(p->min.y - pos.y);
        if (pos.z > p->max.z) mind += sqr(pos.z - p->max.z);
        if (pos.z < p->min.z) mind += sqr(p->min.z - pos.z);
        if (mind > p->maxr2) return;
        if (p->hitpoint->valid && (pos - p->hitpoint->pos).len2() <= p->hitpoint->r2) {
            HitPoint *hp = p->hitpoint;
            double factor = (hp->n * alpha + alpha) / (hp->n * alpha + 1.);
            hp->n++;
            hp->r2 *= factor;
            hp->flux = (hp->flux + hp->fl.mult(fl)) * factor;
        }
        if (p->ls) update(p->ls, pos, fl, d);
        if (p->rs) update(p->rs, pos, fl, d);
        p->maxr2 = p->hitpoint->r2;
        if (p->ls && p->ls->maxr2 > p->maxr2)
            p->maxr2 = p->ls->maxr2;
        if (p->rs && p->rs->maxr2 > p->maxr2)
            p->maxr2 = p->rs->maxr2;
    }
};

HitPointKDTree* hitpointsKDTree;

void initializeHitpointKDTree(std::vector<HitPoint*>& hitpoints) {
    if (hitpointsKDTree)
        delete hitpointsKDTree;
    hitpointsKDTree = new HitPointKDTree(hitpoints);
    // fprintf(stderr, "Hitpoint KD tree built\n");
}

#endif // __HITPOINT_H__