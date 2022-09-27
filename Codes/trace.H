#ifndef TRACE_H
#define TRACE_H

#include "slVector.H"
#include <vector>

class Ray {
public:
    SlVector3 e;
    SlVector3 d;
    int depth;
    Ray(const SlVector3 &_e, const SlVector3 &_d, int _depth = 0) : e(_e), d(_d), depth(_depth) {};
};
class Fill {
public: 
    SlVector3 color;
    double kd, ks, shine, t, ior;
};

class HitRecord {
public:
    double t, alpha, beta, gamma;
    SlVector3 p, n, v;
    Fill f;
    int raydepth;
};

class Light {
public:
    SlVector3 p, c;
};

class Surface {
public:
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const = 0;
    virtual SlVector3 Center() = 0;
    virtual SlVector3 minP() = 0;
    virtual SlVector3 maxP() = 0;
    virtual ~Surface() {};
};

class Triangle : public Surface {
    SlVector3 a,b,c;
public:
    Triangle(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c) : a(_a), b(_b), c(_c) {};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    SlVector3 Center();
    SlVector3 minP();
    SlVector3 maxP();
};

class TrianglePatch : public Triangle {
    SlVector3 n1, n2, n3;
public: 
    TrianglePatch(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c,
    const SlVector3 &_n1, const SlVector3 &_n2, const SlVector3 &_n3) 
        : Triangle(_a,_b,_c), n1(_n1), n2(_n2), n3(_n3) {};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
};

class Sphere : public Surface {
    SlVector3 c;
    double rad;
public:
    Sphere(const SlVector3 &_c, double _r) : c(_c), rad(_r) {};
    bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    SlVector3 Center();
    SlVector3 minP();
    SlVector3 maxP();
};

class AABB {
public:
    AABB() = default;
    bool intersect(const Ray& r, double t0, double t1, HitRecord& hr) const;
    SlVector3 MinC;
    SlVector3 MaxC;
};
struct Node {
    Node() = default;
    Node(Node &node) {
        Box = node.Box;
        lNode = node.lNode;
        rNode = node.rNode;
        PushSur = node.PushSur;
    }
    AABB* Box;
    Node* lNode = nullptr;
    Node* rNode = nullptr;
    std::vector<std::pair<Surface*, Fill> > PushSur;
};

Node* StructNode(const std::vector<std::pair<Surface*, Fill> >& surfaces);

void FreeNode(Node* node) {
    if (!node) return;
    free(node->Box);
    FreeNode(node->lNode);
    FreeNode(node->rNode);
}
class Tracer {
    SlVector3 bcolor, eye, at, up;
    double angle, hither;
    std::vector<Light> lights;
    double shadowbias;
    unsigned int res[2];
  
    SlVector3 *im;
public:
    Tracer(const std::string &fname);
    ~Tracer();
    void traceImage();  
    SlVector3 trace(const Ray &ray, double t0, double t1, const int& _depth) const;

    SlVector3 shade(const HitRecord &hr, const int& _depth) const;

    void writeImage(const std::string &fname);
    bool BVHIntersect(Node* node, const Ray& r, double t0, double t1, HitRecord& hr, bool isShadow = false) const;
    void BVHTransfer();
    double Fresnel_Refract(const SlVector3& I, const SlVector3& N, const double& ior) const;
    SlVector3 Refract(const SlVector3& I, const SlVector3& N, const double& ior) const;

    bool color;
    int samples;
    double aperture;
    int maxraydepth;
    Node* BVH = nullptr;
    std::vector<std::pair<Surface*, Fill> > surfaces;
};

#endif
