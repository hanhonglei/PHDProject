#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <set>

using std::cout;
using std::wcout;
using std::ifstream;
using std::getline;
using std::string;
using std::wstring;
using std::stringstream;
using std::vector;

#include "Color.h"
#include "bsphere.h"
#include "lineqn.h"
#include "Vec.h"

class Mesh
{
public:
    Mesh(){}
    ~Mesh(){}

    // Types
    struct Face {
        vector<point>::size_type v[3];

        Face() {}
        Face(const vector<point>::size_type &v0, 
             const vector<point>::size_type &v1,
             const vector<point>::size_type &v2)
        { v[0] = v0; v[1] = v1; v[2] = v2; }

        Face(const vector<point>::size_type *v_)
        { v[0] = v_[0]; v[1] = v_[1]; v[2] = v_[2]; }

        vector<point>::size_type &operator[] (int i) { return v[i]; }
        const vector<point>::size_type &operator[] (int i) const { return v[i]; }
        operator const vector<point>::size_type * () const { return &(v[0]); }
        operator const vector<point>::size_type * () { return &(v[0]); }
        operator vector<point>::size_type * () { return &(v[0]); }
        int indexof(vector<point>::size_type v_) const
        {
            return (v[0] == v_) ? 0 :
                (v[1] == v_) ? 1 :
                (v[2] == v_) ? 2 : -1;
        }
    };

    struct BBox {
        point min, max;
        point center() const { return 0.5f * (min+max); }
        vec size() const { return max - min; }
    };

    struct BSphere {
        point center;
        float r;
        bool valid;
        BSphere() : valid(false)
        { }
    };

    // per-vertex property
    vector<point> vertices;
    vector<vector<point>::size_type> vertices_sequence;
    vector<Face> faces;
    vector<unsigned int> segmentation;

    vector<Color> colors;
    vector<float> confidences;

    // bounding structures
    BBox bbox;
    BSphere bsphere;

    // computed properties
    vector<vec> normals;
    vector<vec> cornerareas;
    vector<float> pointareas;

    // principal curvature per vertices,
    vector<vec> pdir1, pdir2;
    vector<float> curv1, curv2;
    vector< Vec<4,float> > dcurv;

    // algorithms
    static bool read_helper(const wstring filename, Mesh *mesh);
    static Mesh* read(const wstring filename);
    bool read_segmentation(const wstring filename);

    void compute_bsphere();
    void compute_faces();
    void compute_normals();
    void compute_curvatures();
    void compute_dcurv();
    void compute_pointareas();

    // compute mesh saliency
    std::map<unsigned int, std::set<unsigned int> > adjacent;
    std::vector<float> saliency;
    std::vector<float> saliency_mapped_color; // 3*saliency.size()
    float epsilon; // for considering neighborhood

    float GaussianWeigtedAverage(unsigned int i, float delta);
    void compute_adjacent_vertices();
    void compute_mesh_saliency();
};
