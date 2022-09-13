#include "D:\James O'Brien\tracefiles\trace\trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#define MAX DBL_MAX
#endif
#define M_PI 3.14159265
// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3& a, const SlVector3& b, const SlVector3& c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
        b[0] * (c[1] * a[2] - a[1] * c[2]) +
        c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) { return x * x; }

double Tracer::Fresnel_Refract(const SlVector3& I, const SlVector3& N, const double& ior) const {
    double CosTheta = dot(I, N);
    double n1 = 1.0, n2 = ior;
    double kr;
    if (CosTheta < 0) std::swap(n1, n2);
    double SinTheta = n1 / n2 * std::sqrt(1.0 - sqr(CosTheta));
    if (SinTheta >= 1.0) kr = 1.0;
    else {
        double CosThetat = std::sqrt(1.0 - sqr(SinTheta));
        CosTheta = fabs(CosTheta);
        double FrPara = std::pow((n2 * CosTheta - n1 * CosThetat) / (n2 * CosTheta + n1 * CosThetat), 2.0);
        double FrPerp = std::pow((n1 * CosThetat - n2 * CosTheta) / (n1 * CosThetat + n2 * CosTheta), 2.0);
        kr = (FrPara + FrPerp) / 2.0;
    }
    return kr;
}

SlVector3 Tracer::Refract(const SlVector3& I, const SlVector3& N, const double& ior) const {
    double CosTheta = dot(I, N);
    SlVector3 normal = N;
    double n1 = 1.0, n2 = ior;

    if (CosTheta < 0) {
        std::swap(n1, n2);
        normal = -normal;
        CosTheta = -CosTheta;
    }

    double n1_over_n2 = n1 / n2;
    double discriminant = 1.0 - sqr(n1_over_n2) * (1.0 - sqr(CosTheta));

    if (discriminant < 0) {
        return (0, 0, 0);
    }
    else {
        return - n1_over_n2 * I + (n1_over_n2 * CosTheta - std::sqrt(discriminant)) * normal;
    }
}

bool Triangle::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {

    // Step 1 Ray-triangle test
    double matrix[3][3], Tmatrix[3][3];
    SlVector3 V2subV1 = this->b - this->a;
    SlVector3 V3subV1 = this->c - this->a;
    SlVector3 AsubV1 = r.e - this->a;

    // Original Code with Violent Matrix Calculation Replaced by Prof Code
    /*matrix[0][0] = V2subV1[0];
    matrix[0][1] = V2subV1[1];
    matrix[0][2] = V2subV1[2];
    matrix[1][0] = V3subV1[0];
    matrix[1][1] = V3subV1[1];
    matrix[1][2] = V3subV1[2];
    matrix[2][0] = r.d[0];
    matrix[2][1] = r.d[1];
    matrix[2][2] = r.d[2];
    double Determinant = det(-V2subV1, -V3subV1, r.d);
    if (Determinant == 0) return false;
    Tmatrix[0][0] = (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]) * (1.0 / Determinant);
    Tmatrix[0][1] = (- matrix[1][0] * matrix[2][2] + matrix[2][0] * matrix[1][2]) * (1.0 / Determinant);
    Tmatrix[0][2] = (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]) * (1.0 / Determinant);
    Tmatrix[1][0] = (- matrix[0][1] * matrix[2][2] + matrix[2][1] * matrix[0][2]) * (1.0 / Determinant);
    Tmatrix[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * (1.0 / Determinant);
    Tmatrix[1][2] = (- matrix[0][0] * matrix[2][1] + matrix[2][0] * matrix[0][1]) * (1.0 / Determinant);
    Tmatrix[2][0] = (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) * (1.0 / Determinant);
    Tmatrix[2][1] = (- matrix[0][0] * matrix[1][2] + matrix[0][2] * matrix[1][0]) * (1.0 / Determinant);
    Tmatrix[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * (1.0 / Determinant);
    double Beta = Tmatrix[0][0] * AsubV1[0] + Tmatrix[0][1] * AsubV1[1] + Tmatrix[0][2] * AsubV1[2];
    double Gamma = Tmatrix[1][0] * AsubV1[0] + Tmatrix[1][1] * AsubV1[1] + Tmatrix[1][2] * AsubV1[2];
    double T = -(Tmatrix[2][0] * AsubV1[0] + Tmatrix[2][1] * AsubV1[1] + Tmatrix[2][2] * AsubV1[2]);*/

    double Determinant = det(-V2subV1, -V3subV1, r.d);
    double T = det(-V2subV1, -V3subV1, -AsubV1) / Determinant;
    double Beta = det(-AsubV1, -V3subV1, r.d) / Determinant;
    double Gamma = det(-V2subV1, -AsubV1, r.d) / Determinant;

    if (T < t0 || T > t1) return false;
    if (Beta < 0.0 || Beta > 1.0) return false;
    else if (Gamma < 0.0 || Gamma > 1.0 - Beta) return false;
    else {
        hr.beta = Beta;
        hr.gamma = Gamma;
        hr.t = T;
        hr.p = r.e + T * r.d;
        hr.alpha = 1.0 - Beta - Gamma;
        hr.n = cross(V2subV1, V3subV1);
        normalize(hr.n);
        return true;
    }
}

bool TrianglePatch::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {
    bool temp = Triangle::intersect(r, t0, t1, hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
    }
    return temp;
}

bool Sphere::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {

    // Step 1 Sphere-triangle test
    double D_AC = dot(r.d, r.e - this->c);
    double D_Square = sqrMag(r.d);
    double AsubC = sqrMag(r.e - this->c);
    double delta = sqr(D_AC) - D_Square * (AsubC - sqr(this->rad));
    if (delta < 0.0) return false;
    else {
        double t_1 = (- D_AC + sqrt(delta)) / D_Square;
        double t_2 = (- D_AC - sqrt(delta)) / D_Square;
        double tmin = t_1;
        if (t_1 < 0 || (t_2 > 0 && t_1 > t_2)) tmin = t_2;
        if (tmin < t0 || tmin > t1) return false;
        hr.t = tmin;
        hr.p = r.e + tmin * r.d;
        hr.n = (hr.p - this->c) / this->rad;
        return true;
    }
}

bool AABB::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {
    this->MinC;
    this->MaxC;
    double tmin = (MinC.x() - r.e.x()) / r.d.x();
    double tmax = (MaxC.x() - r.e.x()) / r.d.x();

    if (tmin > tmax) std::swap(tmin, tmax);

    double tymin = (MinC.y() - r.e.y()) / r.d.y();
    double tymax = (MaxC.y() - r.e.y()) / r.d.y();

    if (tymin > tymax) std::swap(tymin, tymax);

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    tmin = std::min(tmin, tymin);
    tmax = std::max(tmax, tymax);

    double tzmin = (MinC.z() - r.e.z()) / r.d.z();
    double tzmax = (MaxC.z() - r.e.z()) / r.d.z();

    if (tzmin > tzmax) std::swap(tzmin, tzmax);

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    tmin = std::min(tmin, tzmin);
    tmax = std::max(tmax, tzmax);

    if (tmax < t0 || tmin > t1) return false;
    hr.t = tmin;

    return true;
}

SlVector3 Sphere::Center() {
    return c;
}
SlVector3 Triangle::Center() {
    return (a + b + c) / 3;
}
SlVector3 Sphere::minP() {
    return {
        c.x() - rad,
        c.y() - rad,
        c.z() - rad
    };
}
SlVector3 Sphere::maxP() {
    return {
        c.x() + rad,
        c.y() + rad,
        c.z() + rad
    };
}
SlVector3 Triangle::minP() {
    return {
        std::min(std::min(a.x(), b.x()), c.x()),
        std::min(std::min(a.y(), b.y()), c.y()),
        std::min(std::min(a.z(), b.z()), c.z())
    };
}
SlVector3 Triangle::maxP() {
    return {
        std::max(std::max(a.x(), b.x()), c.x()),
        std::max(std::max(a.y(), b.y()), c.y()),
        std::max(std::max(a.z(), b.z()), c.z())
    };
}

Node* StructNode(const std::vector<std::pair<Surface *, Fill> >& surfaces) {
    if (surfaces.empty()) return nullptr;

    Node* node = new Node();
    node->PushSur = surfaces;

    AABB* Box = new AABB();
    Box->MinC = surfaces[0].first->minP();
    Box->MaxC = surfaces[0].first->maxP();
    SlVector3 MinCenter, MaxCenter;
    MinCenter = MaxCenter = surfaces[0].first->Center();

    for (unsigned int i = 1; i < surfaces.size(); i++) {
        Surface* surface = surfaces[i].first;
        MinCenter = min(MinCenter, surface->Center());
        MaxCenter = max(MaxCenter, surface->Center());
        Box->MinC = min(Box->MinC, surface->minP());
        Box->MaxC = max(Box->MaxC, surface->maxP());
    }
    node->Box = Box;
    if (surfaces.size() == 1) return node;

    std::vector<std::pair<Surface*, Fill> > leftNode, rightNode;
    SlVector3 difference = MaxCenter - MinCenter;
    unsigned int Index = 0;
    for (int i = 1; i < 3; i++) {
        if (difference[i] > difference[Index]) Index = i;
    }
    double threshold = (MinCenter[Index] + MaxCenter[Index]) / 2.0;
    for (unsigned int i = 0; i < surfaces.size(); i ++) {
        const std::pair<Surface*, Fill>& s = surfaces[i];
        if (s.first->Center()[Index] < threshold) leftNode.push_back(s);
        else rightNode.push_back(s);
    }

    node->lNode = StructNode(leftNode);
    node->rNode = StructNode(rightNode);

    return node;
}

Tracer::Tracer(const std::string& fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        if (in.eof()) break;
        getline(in, line);
        switch (line[0]) {
        case 'b': {
            std::stringstream ss(line);
            ss >> ch >> bcolor[0] >> bcolor[1] >> bcolor[2];
            break;
        }

        case 'v': {
            getline(in, line);
            std::string junk;
            std::stringstream fromss(line);
            fromss >> junk >> eye[0] >> eye[1] >> eye[2];

            getline(in, line);
            std::stringstream atss(line);
            atss >> junk >> at[0] >> at[1] >> at[2];

            getline(in, line);
            std::stringstream upss(line);
            upss >> junk >> up[0] >> up[1] >> up[2];

            getline(in, line);
            std::stringstream angless(line);
            angless >> junk >> angle;

            getline(in, line);
            std::stringstream hitherss(line);
            hitherss >> junk >> hither;

            getline(in, line);
            std::stringstream resolutionss(line);
            resolutionss >> junk >> res[0] >> res[1];
            break;
        }

        case 'p': {
            bool patch = false;
            std::stringstream ssn(line);
            unsigned int nverts;
            if (line[1] == 'p') {
                patch = true;
                ssn >> ch;
            }
            ssn >> ch >> nverts;
            std::vector<SlVector3> vertices;
            std::vector<SlVector3> normals;
            for (unsigned int i = 0; i < nverts; i++) {
                getline(in, line);
                std::stringstream ss(line);
                SlVector3 v, n;
                if (patch) ss >> v[0] >> v[1] >> v[2] >> n[0] >> n[1] >> n[2];
                else ss >> v[0] >> v[1] >> v[2];
                vertices.push_back(v);
                normals.push_back(n);
            }
            bool makeTriangles = false;
            if (vertices.size() == 3) {
                if (patch) {
                    surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                        normals[0], normals[1], normals[2]), fill));
                }
                else {
                    surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                }
            }
            else if (vertices.size() == 4) {
                SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                if (dot(n0, n1) > 0 && dot(n0, n2) > 0 && dot(n0, n3) > 0) {
                    makeTriangles = true;
                    if (patch) {
                        surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                            normals[0], normals[1], normals[2]), fill));
                        surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3],
                            normals[0], normals[2], normals[3]), fill));
                    }
                    else {
                        surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                        surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
                    }
                }
                if (!makeTriangles) {
                    std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                }
            }
            break;
        }

        case 's': {
            std::stringstream ss(line);
            SlVector3 c;
            double r;
            ss >> ch >> c[0] >> c[1] >> c[2] >> r;
            surfaces.push_back(std::pair<Surface*, Fill>(new Sphere(c, r), fill));
            break;
        }

        case 'f': {
            std::stringstream ss(line);
            ss >> ch >> fill.color[0] >> fill.color[1] >> fill.color[2] >> fill.kd >> fill.ks >> fill.shine >> fill.t >> fill.ior;
            break;
        }

        case 'l': {
            std::stringstream ss(line);
            Light l;
            ss >> ch >> l.p[0] >> l.p[1] >> l.p[2];
            if (!ss.eof()) {
                ss >> l.c[0] >> l.c[1] >> l.c[2];
                coloredlights = true;
            }
            lights.push_back(l);
            break;
        }

        default:
            break;
        }
    }
    if (!coloredlights) for (unsigned int i = 0; i < lights.size(); i++) lights[i].c = 1.0 / sqrt(lights.size());
    im = new SlVector3[res[0] * res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}

Tracer::~Tracer() {
    if (im) delete[] im;
    for (unsigned int i = 0; i < surfaces.size(); i++) delete surfaces[i].first;
    FreeNode(this->BVH);
}

SlVector3 Tracer::shade(const HitRecord& hr, const int& _depth) const {
    if (color) return hr.f.color;

    SlVector3 color(0.0), RefrColor(0.0), ReflColor(0.0);
    HitRecord dummy;
    for (unsigned int i = 0; i < lights.size(); i++) {
        const Light& light = lights[i];
        bool shadow = false;

        // Step 3 Check for shadows here
        SlVector3 SurLight = light.p - hr.p;
        Ray SurRay (hr.p, SurLight);
        normalize(SurRay.d);
        shadow = BVHIntersect(BVH, SurRay, this->shadowbias, mag(SurLight), dummy, true);

        if (!shadow) {
            // Step 2 do shading here
            SlVector3 LightRay = light.p - hr.p;
            normalize (LightRay);
            SlVector3 ReflectionRay = -LightRay + 2 * dot(LightRay, hr.n) * hr.n;
            normalize (ReflectionRay);
            SlVector3 DiffuseColor, SpecularColor, AmbientColor;
            DiffuseColor = light.c * hr.f.kd * std::max(dot(LightRay, hr.n), 0.0) * hr.f.color;
            SpecularColor = light.c * hr.f.ks * std::pow(std::max(dot(hr.v, ReflectionRay), 0.0), hr.f.shine) * hr.f.color;
            AmbientColor = light.c * 0.01 * hr.f.color;
            color += DiffuseColor + SpecularColor + AmbientColor;
        }
    }

    // Step 4 Add code for computing reflection color here
    if (_depth > this->maxraydepth) return{};

    if (hr.f.ks > 0) {
        SlVector3 DirRay = -hr.v + 2 * dot(hr.v, hr.n) * hr.n;
        Ray Refray(hr.p, DirRay);
        normalize(Refray.d);
        color += .5 * (hr.f.ks + hr.f.kd) * hr.f.ks * trace(Refray, this->shadowbias, MAX, _depth + 1);
    }

    // Step 5 Add code for computing refraction color here
    if (hr.f.t > 0) {
        double kr = this->Fresnel_Refract(hr.v, hr.n, hr.f.ior);
        //SlVector3 RefractColor = (0, 0, 0);
        SlVector3 RefractDir = this->Refract(hr.v, hr.n, hr.f.ior);

        if (kr < 1) {
            Ray RefractRay(hr.p, RefractDir);
            normalize(RefractRay.d);
            RefrColor += trace(RefractRay, this->shadowbias, MAX, _depth + 1) * hr.f.color;
        }

        SlVector3 normal = hr.n;
        if (RefractDir == (0, 0, 0)) {
            RefractDir = hr.v;
            normal = -normal;
        }
        SlVector3 ReflectDir = -RefractDir + 2 * dot(RefractDir, normal) * normal;
        Ray ReflectRay(hr.p, ReflectDir);
        //SlVector3 ReflectColor = trace(ReflectRay, this->shadowbias, MAX, _depth + 1);
        ReflColor += trace(ReflectRay, this->shadowbias, MAX, _depth + 1) * hr.f.color;
        color += ReflColor * kr + RefrColor * (1.0 - kr);
    }
    
    return color;
}

SlVector3 Tracer::trace(const Ray& r, double t0, double t1, const int& _depth) const {
    HitRecord hr;
    SlVector3 color(bcolor);
    bool hit = false;

    // Step 1 See what a ray hits 
    hit = BVHIntersect(BVH, r, t0, t1, hr);
    
    if (hit) color = shade(hr, _depth);
    
    return color;
}

bool Tracer::BVHIntersect(Node* node, const Ray& r, double t0, double t1, HitRecord& hr, bool isShadow) const {
    bool hit = false;
    HitRecord Hit;
    std::stack<Node*> VStack;
    VStack.push(node);

    while (!VStack.empty()) {
        if (isShadow && hit) return true;
        Node* CurrentNode = VStack.top();
        VStack.pop();
        if (CurrentNode->Box->intersect(r, t0, t1, Hit)) {
            if (CurrentNode->PushSur.size() == 1) {
                const std::pair<Surface*, Fill>& s = CurrentNode->PushSur[0];
                if (s.first->intersect(r, t0, t1, hr)) {
                    hit = true;
                    t1 = hr.t;
                    hr.raydepth = r.depth;
                    hr.v = r.e - hr.p;
                    hr.f = s.second;
                    normalize(hr.v);
                    normalize(hr.n);
                }
            }
            if (CurrentNode->lNode) VStack.push(CurrentNode->lNode);
            if (CurrentNode->rNode) VStack.push(CurrentNode->rNode);
        }
    }
    return hit;
}

void Tracer::traceImage() {
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up, w);
    normalize(u);
    SlVector3 v = cross(w, u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI / 180.0) * (angle / 2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3* pixel = im;

    for (unsigned int j = 0; j < res[1]; j++) {
        for (unsigned int i = 0; i < res[0]; i++, pixel++) {

            SlVector3 result(0.0, 0.0, 0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r - l) * (i + rx) / res[0];
                double y = b + (t - b) * (j + ry) / res[1];
                SlVector3 dir = -d * w + x * u + y * v;

                Ray r(eye, dir);
                normalize(r.d);

                result += trace(r, hither, MAX, 0);

            }
            (*pixel) = result / samples;
        }
    }
}

void Tracer::writeImage(const std::string& fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out << "P6" << "\n" << res[0] << " " << res[1] << "\n" << 255 << "\n";
    SlVector3* pixel = im;
    char val;
    for (unsigned int i = 0; i < res[0] * res[1]; i++, pixel++) {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write(&val, sizeof(unsigned char));
    }
    out.close();
}

void Tracer::BVHTransfer() {
    BVH = StructNode(surfaces);
}

int main(int argc, char* argv[]) {
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    Tracer tracer("D:\\James O'Brien\\tracefiles\\as2-InputFiles\\InputFiles\\refract.nff");
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    tracer.BVHTransfer();
    tracer.traceImage();
    tracer.writeImage("D:\\James O'Brien\\tracefiles\\as2-InputFiles\\InputFiles\\refract.ppm");
};
