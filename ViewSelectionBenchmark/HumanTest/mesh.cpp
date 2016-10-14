#include "mesh.h"

void Mesh::compute_bsphere()
{
    if (vertices.empty() || bsphere.valid)
        return;
    Miniball<3,float> mb;
    mb.check_in(vertices.begin(), vertices.end());
    mb.build();
    bsphere.center = mb.center();
    bsphere.r = sqrt(mb.squared_radius());
    bsphere.valid = true; 
}

// i+1 and i-1 modulo 3
// This way of computing it tends to be faster than using %
#define NEXT(i) ((i)<2 ? (i)+1 : (i)-2)
#define PREV(i) ((i)>0 ? (i)-1 : (i)+2)

// Rotate a coordinate system to be perpendicular to the given normal
static void rot_coord_sys(const vec &old_u, const vec &old_v,
                          const vec &new_norm,
                          vec &new_u, vec &new_v)
{
    new_u = old_u;
    new_v = old_v;
    vec old_norm = old_u CROSS old_v;
    float ndot = old_norm DOT new_norm;
    if (unlikely(ndot <= -1.0f)) {
        new_u = -new_u;
        new_v = -new_v;
        return;
    }
    vec perp_old = new_norm - ndot * old_norm;
    vec dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
    new_u -= dperp * (new_u DOT perp_old);
    new_v -= dperp * (new_v DOT perp_old);
}

// Reproject a curvature tensor from the basis spanned by old_u and old_v
// (which are assumed to be unit-length and perpendicular) to the
// new_u, new_v basis.
static void proj_curv(const vec &old_u, const vec &old_v,
               float old_ku, float old_kuv, float old_kv,
               const vec &new_u, const vec &new_v,
               float &new_ku, float &new_kuv, float &new_kv)
{
    vec r_new_u, r_new_v;
    rot_coord_sys(new_u, new_v, old_u CROSS old_v, r_new_u, r_new_v);

    float u1 = r_new_u DOT old_u;
    float v1 = r_new_u DOT old_v;
    float u2 = r_new_v DOT old_u;
    float v2 = r_new_v DOT old_v;
    new_ku  = old_ku * u1*u1 + old_kuv * (2.0f  * u1*v1) + old_kv * v1*v1;
    new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
    new_kv  = old_ku * u2*u2 + old_kuv * (2.0f  * u2*v2) + old_kv * v2*v2;
}

// Like the above, but for dcurv
static void proj_dcurv(const vec &old_u, const vec &old_v,
                const Vec<4> old_dcurv,
                const vec &new_u, const vec &new_v,
                Vec<4> &new_dcurv)
{
    vec r_new_u, r_new_v;
    rot_coord_sys(new_u, new_v, old_u CROSS old_v, r_new_u, r_new_v);

    float u1 = r_new_u DOT old_u;
    float v1 = r_new_u DOT old_v;
    float u2 = r_new_v DOT old_u;
    float v2 = r_new_v DOT old_v;

    new_dcurv[0] = old_dcurv[0]*u1*u1*u1 +
        old_dcurv[1]*3.0f*u1*u1*v1 +
        old_dcurv[2]*3.0f*u1*v1*v1 +
        old_dcurv[3]*v1*v1*v1;
    new_dcurv[1] = old_dcurv[0]*u1*u1*u2 +
        old_dcurv[1]*(u1*u1*v2 + 2.0f*u2*u1*v1) +
        old_dcurv[2]*(u2*v1*v1 + 2.0f*u1*v1*v2) +
        old_dcurv[3]*v1*v1*v2;
    new_dcurv[2] = old_dcurv[0]*u1*u2*u2 +
        old_dcurv[1]*(u2*u2*v1 + 2.0f*u1*u2*v2) +
        old_dcurv[2]*(u1*v2*v2 + 2.0f*u2*v2*v1) +
        old_dcurv[3]*v1*v2*v2;
    new_dcurv[3] = old_dcurv[0]*u2*u2*u2 +
        old_dcurv[1]*3.0f*u2*u2*v2 +
        old_dcurv[2]*3.0f*u2*v2*v2 +
        old_dcurv[3]*v2*v2*v2;
}


// Given a curvature tensor, find principal directions and curvatures
// Makes sure that pdir1 and pdir2 are perpendicular to normal
static void diagonalize_curv(const vec &old_u, const vec &old_v,
                      float ku, float kuv, float kv,
                      const vec &new_norm,
                      vec &pdir1, vec &pdir2, float &k1, float &k2)
{
    vec r_old_u, r_old_v;
    rot_coord_sys(old_u, old_v, new_norm, r_old_u, r_old_v);

    float c = 1, s = 0, tt = 0;
    if (likely(kuv != 0.0f)) {
        // Jacobi rotation to diagonalize
        float h = 0.5f * (kv - ku) / kuv;
        tt = (h < 0.0f) ?
            1.0f / (h - sqrt(1.0f + h*h)) :
        1.0f / (h + sqrt(1.0f + h*h));
        c = 1.0f / sqrt(1.0f + tt*tt);
        s = tt * c;
    }

    k1 = ku - tt * kuv;
    k2 = kv + tt * kuv;

    if (fabs(k1) >= fabs(k2)) {
        pdir1 = c*r_old_u - s*r_old_v;
    } else {
        swap(k1, k2);
        pdir1 = s*r_old_u + c*r_old_v;
    }
    pdir2 = new_norm CROSS pdir1;
}

// Compute principal curvatures and directions.
void Mesh::compute_curvatures()
{
    if (curv1.size() == vertices.size())
        return;
    compute_faces();
    compute_normals();
    compute_pointareas();

    // Resize the arrays we'll be using
    vector<Face>::size_type nf = faces.size();
    vector<point>::size_type nv = vertices.size();

    curv1.clear(); curv1.resize(nv); curv2.clear(); curv2.resize(nv);
    pdir1.clear(); pdir1.resize(nv); pdir2.clear(); pdir2.resize(nv);
    vector<float> curv12(nv);

    // Set up an initial coordinate system per vertex
    for (vector<Face>::size_type i = 0; i < nf; i++) {
        pdir1[faces[i][0]] = vertices[faces[i][1]] -
            vertices[faces[i][0]];
        pdir1[faces[i][1]] = vertices[faces[i][2]] -
            vertices[faces[i][1]];
        pdir1[faces[i][2]] = vertices[faces[i][0]] -
            vertices[faces[i][2]];
    }
    for (vector<point>::size_type i = 0; i < nv; i++) {
        pdir1[i] = pdir1[i] CROSS normals[i];
        normalize(pdir1[i]);
        pdir2[i] = normals[i] CROSS pdir1[i];
    }

    // Compute curvature per-face
    for (vector<Face>::size_type i = 0; i < nf; i++) {
        // Edges
        vec e[3] = { vertices[faces[i][2]] - vertices[faces[i][1]],
            vertices[faces[i][0]] - vertices[faces[i][2]],
            vertices[faces[i][1]] - vertices[faces[i][0]] };

        // N-T-B coordinate system per face
        vec t = e[0];
        normalize(t);
        vec n = e[0] CROSS e[1];
        vec b = n CROSS t;
        normalize(b);

        // Estimate curvature based on variation of normals
        // along edges
        float m[3] = { 0, 0, 0 };
        float w[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };
        for (int j = 0; j < 3; j++) {
            float u = e[j] DOT t;
            float v = e[j] DOT b;
            w[0][0] += u*u;
            w[0][1] += u*v;
            //w[1][1] += v*v + u*u; 
            //w[1][2] += u*v; 
            w[2][2] += v*v;
            vec dn = normals[faces[i][PREV(j)]] -
                normals[faces[i][NEXT(j)]];
            float dnu = dn DOT t;
            float dnv = dn DOT b;
            m[0] += dnu*u;
            m[1] += dnu*v + dnv*u;
            m[2] += dnv*v;
        }
        w[1][1] = w[0][0] + w[2][2];
        w[1][2] = w[0][1];

        // Least squares solution
        float diag[3];
        if (!ldltdc<float,3>(w, diag)) {
            //fprintf(stderr, "ldltdc failed!\n");
            continue;
        }
        ldltsl<float,3>(w, diag, m, m);

        // Push it back out to the vertices
        for (int j = 0; j < 3; j++) {
            vector<point>::size_type vj = faces[i][j];
            float c1, c12, c2;
            proj_curv(t, b, m[0], m[1], m[2],
                pdir1[vj], pdir2[vj], c1, c12, c2);
            float wt = cornerareas[i][j] / pointareas[vj];
            curv1[vj]  += wt * c1;
            curv12[vj] += wt * c12;
            curv2[vj]  += wt * c2;
        }
    }

    // Compute principal directions and curvatures at each vertex
    for (vector<point>::size_type i = 0; i < nv; i++)
        diagonalize_curv(pdir1[i], pdir2[i],
        curv1[i], curv12[i], curv2[i],
        normals[i], pdir1[i], pdir2[i],
        curv1[i], curv2[i]);
}

// Compute derivatives of curvature.
void Mesh::compute_dcurv()
{
    if (dcurv.size() == vertices.size())
        return;
    compute_curvatures();

    // Resize the arrays we'll be using
    vector<Face>::size_type nf = faces.size();
    vector<point>::size_type nv = vertices.size();
    dcurv.clear(); dcurv.resize(nv);

    // Compute dcurv per-face
    for (vector<Face>::size_type i = 0; i < nf; i++) {
        // Edges
        vec e[3] = { vertices[faces[i][2]] - vertices[faces[i][1]],
            vertices[faces[i][0]] - vertices[faces[i][2]],
            vertices[faces[i][1]] - vertices[faces[i][0]] };

        // N-T-B coordinate system per face
        vec t = e[0];
        normalize(t);
        vec n = e[0] CROSS e[1];
        vec b = n CROSS t;
        normalize(b);

        // Project curvature tensor from each vertex into this
        // face's coordinate system
        vec fcurv[3];
        for (int j = 0; j < 3; j++) {
            vector<point>::size_type vj = faces[i][j];
            proj_curv(pdir1[vj], pdir2[vj], curv1[vj], 0, curv2[vj],
                t, b, fcurv[j][0], fcurv[j][1], fcurv[j][2]);

        }

        // Estimate dcurv based on variation of curvature along edges
        float m[4] = { 0, 0, 0, 0 };
        float w[4][4] = { {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0} };
        for (int j = 0; j < 3; j++) {
            // Variation of curvature along each edge
            vec dfcurv = fcurv[PREV(j)] - fcurv[NEXT(j)];
            float u = e[j] DOT t;
            float v = e[j] DOT b;
            float u2 = u*u, v2 = v*v, uv = u*v;
            w[0][0] += u2;
            w[0][1] += uv;
            //w[1][1] += 2.0f*u2 + v2;
            //w[1][2] += 2.0f*uv;
            //w[2][2] += u2 + 2.0f*v2;
            //w[2][3] += uv;
            w[3][3] += v2;
            m[0] += u*dfcurv[0];
            m[1] += v*dfcurv[0] + 2.0f*u*dfcurv[1];
            m[2] += 2.0f*v*dfcurv[1] + u*dfcurv[2];
            m[3] += v*dfcurv[2];
        }
        w[1][1] = 2.0f * w[0][0] + w[3][3];
        w[1][2] = 2.0f * w[0][1];
        w[2][2] = w[0][0] + 2.0f * w[3][3];
        w[2][3] = w[0][1];

        // Least squares solution
        float d[4];
        if (!ldltdc<float,4>(w, d)) {
            //fprintf(stderr, "ldltdc failed!\n");
            continue;
        }
        ldltsl<float,4>(w, d, m, m);
        Vec<4> face_dcurv(m);

        // Push it back out to each vertex
        for (int j = 0; j < 3; j++) {
            vector<point>::size_type vj = faces[i][j];
            Vec<4> this_vert_dcurv;
            proj_dcurv(t, b, face_dcurv,
                pdir1[vj], pdir2[vj], this_vert_dcurv);
            float wt = cornerareas[i][j] / pointareas[vj];
            dcurv[vj] += wt * this_vert_dcurv;
        }
    }
}

#define BIGNUM 1.0e38

// Forward declarations
static bool read_ply(ifstream &file, Mesh *mesh);
static bool read_obj(ifstream &file, Mesh *mesh);
static bool read_off(ifstream &file, Mesh *mesh);
static void check_ind_range(Mesh *mesh);

static bool read_ply(ifstream &file, Mesh *mesh)
{
    // read ply file from the second line,
    // ascii 1.0 format only

    string linebuffer,dump;
    stringstream streambuffer;

    vector<point>::size_type nverts;
    vector<Mesh::Face>::size_type nfaces;
    string property_type,property_name;
    vector<string> vertex_property_sequence;
    vector<string> face_property_sequence;

    getline(file,linebuffer); // read "ply" first line
    getline(file,linebuffer); // second line for file format info

    if (linebuffer.find("format ascii 1.0") == string::npos)
    {
        cout<<"only support ascii 1.0 ply file format\n";
        return false;
    }

    // start to parse file header
    while (getline(file,linebuffer),linebuffer != string("end_header"))
    {
        // skipping comments
        if (linebuffer.find("comment") != string::npos) continue;
        if (linebuffer.find("element vertex") != string::npos)
        {
            streambuffer.clear();
            streambuffer.str(linebuffer);
            streambuffer>>dump>>dump>>nverts;

            while (getline(file,linebuffer),
                linebuffer.find("property") != string::npos)
            {
                streambuffer.clear();
                streambuffer.str(linebuffer);
                streambuffer>>dump>>property_type>>property_name;
                vertex_property_sequence.push_back(property_name);
            }
            // putback the next element type info
            file.putback('\n');
            for (string::size_type i=linebuffer.size(); i>0; i--)
                file.putback(linebuffer[i-1]);
        }
        else if (linebuffer.find("element face") != string::npos)
        {
            streambuffer.clear();
            streambuffer.str(linebuffer);
            streambuffer>>dump>>dump>>nfaces;

            while (getline(file,linebuffer),
                linebuffer.find("property list") != string::npos)
            {
                streambuffer.clear();
                streambuffer.str(linebuffer);
                streambuffer>>dump>>dump;
                while (!streambuffer.eof())
                {
                    streambuffer>>property_name;
                    face_property_sequence.push_back(property_name);
                }
            }

            // putback the next element type info
            file.putback('\n');
            for (string::size_type i=linebuffer.size(); i>0; i--)
                file.putback(linebuffer[i-1]);
        }
        else
        {
            cout<<"Unsupported ply element type!\n";
            return false;
        }
    }

    // now, just read position and color data
    // reading date from file
    point tp;   // position
    Color tc;   // color
    bool has_color=false;
    bool has_confidences=false;
    for (vector<string>::size_type i=0; i<vertex_property_sequence.size(); i++)
    {
        if (vertex_property_sequence[i] == "red"    ||
            vertex_property_sequence[i] == "green"  ||
            vertex_property_sequence[i] == "blue")
        {
            has_color = true;
        }
        if (vertex_property_sequence[i] == "confidences")
        {
            has_confidences = true;
        }
    }
    mesh->vertices.resize(nverts);
    if(has_color) mesh->colors.resize(nverts);
    if(has_confidences) mesh->confidences.resize(nverts);

    for (vector<point>::size_type i=0; i<nverts; i++)
    {
        for (vector<string>::size_type j=0; j<vertex_property_sequence.size(); j++)
        {
            if (vertex_property_sequence[j] == "x")
                file>>tp[0];
            else if(vertex_property_sequence[j] == "y")
                file>>tp[1];
            else if(vertex_property_sequence[j] == "z")
                file>>tp[2];
            else if(vertex_property_sequence[j] == "red")
                file>>tc[0];
            else if(vertex_property_sequence[j] == "green")
                file>>tc[1];
            else if(vertex_property_sequence[j] == "blue")
                file>>tc[2];
            else if(vertex_property_sequence[j] == "confidences")
                file>>mesh->confidences[i];
            else
                file>>dump;
        }
        mesh->vertices[i] = tp;
        if(has_color) mesh->colors[i] = tc;
    }   
    cout<<nverts<<" vertices loaded...\n";

    mesh->faces.resize(nfaces);
    int fver;   //vertex per face, dump value
    for (vector<Mesh::Face>::size_type i=0; i<nfaces; i++)
    {
        if(file>>fver, fver != 3)
        {
            cout<<"Unsupported face type!\n";
            return false;
        }
        file>>mesh->faces[i][0]>>mesh->faces[i][1]>>mesh->faces[i][2];
    }
    cout<<nfaces<<" triangles loaded...\n";
    return true;
}

static bool read_obj(ifstream &file, Mesh *mesh)
{
    // read obj file,
    // just support triangle mesh, 
    // and doesn't handle material... etc.
    // only geometry attributions are handled

    string linebuffer,dump;
    stringstream streambuffer;

    // for reading faces attributions list
    bool has_textrue_coord = false;
    bool has_vertex_normal = false;

    char cdump;     //dump char
    string sdump;   //dump string

    while (!file.eof())
    {
        getline(file,linebuffer);
        // skip comments
        if (linebuffer[0] == '#'||
            linebuffer.size() == 0) 
            continue;

        //init stream buffer
        streambuffer.clear();
        streambuffer.str(linebuffer);
        //parse the buffer
        string::size_type first_spc = linebuffer.find_first_of(' ');
        string first_token = linebuffer.substr(0,first_spc);

        if (first_token == string("g") ||
            first_token == string("mtllib") ||
            first_token == string("usemtl"))
        {
            //doesn't handle group name, mtllib, etc
            continue;
        }
        else if (first_token == string("v"))
        {
            //read vertices
            point v;
            streambuffer>>cdump>>v[0]>>v[1]>>v[2];
            mesh->vertices.push_back(v);
        }
        else if (first_token == string("vt"))
        {
            has_textrue_coord = true;

            //read texture coordinate data
            
            //TODO: handle texture coordinates
            //we don't need coordinates in vs project            
        }
        else if (first_token == string("vn"))
        {
            has_vertex_normal = true;
            //read normal data
            //TODO: handle model normal data
            //we don't need model contained normal data in vs project
        }
        else if (first_token == string("f"))
        {
            //read face data
            streambuffer>>cdump;
            Mesh::Face f;
            streambuffer>>f[0];
            if (has_textrue_coord || has_vertex_normal)
                streambuffer>>sdump;
            streambuffer>>f[1];
            if (has_textrue_coord || has_vertex_normal)
                streambuffer>>sdump;
            streambuffer>>f[2];

            //coz the sequence is started from 1, so, sub 1 to get 0 started sequence
            f[0] -= 1;f[1] -= 1;f[2] -= 1;
            mesh->faces.push_back(f);
        }
        else
        {
            cout<<"Unexpected token: "<<first_token<<"\n";
            cout<<"Please check the file again....\n";
            return false;
        }
    }
    cout<<mesh->vertices.size()<<" vertices loaded\n"
        <<mesh->faces.size()<<" faces loaded\n";
    return true;
}

static bool read_off(ifstream &file, Mesh *mesh)
{
    // read off file,
    // just support triangle mesh, 
    // and doesn't handle material... etc.
    // only geometry attributions are handled

    string linebuffer,dump;
    getline(file,linebuffer);
    if (linebuffer != "OFF" && linebuffer != "off")
    {
        std::cout<<"Not off file format...\n";
        return false;
    }

    unsigned int nv = 0;
    unsigned int fv = 0;
    file>>nv>>fv>>dump;

    // read vertices
    mesh->vertices.resize(nv);
    for (unsigned int i=0; i<nv; i++)
    {
        file
            >>mesh->vertices[i][0]
            >>mesh->vertices[i][1]
            >>mesh->vertices[i][2];
    }

    // read faces
    mesh->faces.resize(fv);
    for (unsigned int i=0; i<fv; i++)
    {
        file
            >>dump
            >>mesh->faces[i][0]
            >>mesh->faces[i][1]
            >>mesh->faces[i][2];
    }
    
    cout<<mesh->vertices.size()<<" vertices loaded\n"
        <<mesh->faces.size()<<" faces loaded\n";
    return true;
}

// Check whether the indices in the file mistakenly go
// from 1..N instead of 0..N-1
static void check_ind_range(Mesh *mesh)
{
    if (mesh->faces.empty())
        return;
    vector<point>::size_type min_ind = mesh->faces[0][0];
    vector<point>::size_type max_ind = mesh->faces[0][0];

    for (vector<Mesh::Face>::size_type i = 0; i < mesh->faces.size(); i++) {
        for (int j = 0; j < 3; j++) {
            min_ind = min_ind<mesh->faces[i][j]?min_ind:mesh->faces[i][j];
            max_ind = max_ind>mesh->faces[i][j]?min_ind:mesh->faces[i][j];
        }
    }

    vector<point>::size_type nv = mesh->vertices.size();

    // All good
    if (min_ind == 0 && max_ind == nv-1)
        return;

    // Simple fix: offset everything
    if (max_ind - min_ind == nv-1) {
        cout<<"Found indices ranging from "<<min_ind<<" through "<<max_ind<<"\n";
        cout<<"Remapping to "<<0<<" through "<<nv-1<<"\n";
        for (vector<Mesh::Face>::size_type i = 0; i < mesh->faces.size(); i++)
            for (int j = 0; j < 3; j++)
                mesh->faces[i][j] -= min_ind;
        return;
    }
    // Else can't do anything...
}

Mesh* Mesh::read(const wstring filename)
{
    Mesh *mesh = new Mesh();
    if (read_helper(filename, mesh))
        return mesh;
    delete mesh;
    return NULL;
}

// Actually read a mesh.  Tries to figure out type of file from first
// few bytes.  Filename can be "-" for stdin
bool Mesh::read_helper(const wstring filename, Mesh *mesh)
{
    if (filename.empty())
        return false;

    ifstream file;
    bool ok = false;
    wstring ext(filename);   //file extension
    size_t iext = ext.find_last_of(L'.');
    ext = ext.substr(iext+1,ext.size());

    file.open(filename.c_str());
    if (file.fail()) {
        cout<<"file open error in Mesh:read_helper()\n";
        goto out;
    }

    if (ext == wstring(L"ply"))
    {
        ok = read_ply(file,mesh);
    }
    else if (ext == wstring(L"obj"))
    {
        ok = read_obj(file,mesh);
    }
    else if (ext == wstring(L"off"))
    {
        ok = read_off(file,mesh);
    }
    else
    {
        cout<<"Unknown file type\n";
    }

out:
    file.close();
    if (!ok || (mesh->vertices.empty() && mesh->faces.empty())) 
    {
        cout<<"\nError in reading file ...\n";
        return false;
    }
    check_ind_range(mesh);
    return true;
}

void Mesh::compute_normals()
{
    if (normals.size() == vertices.size())
        return;
    compute_faces();

    normals.clear();
    normals.resize(vertices.size());

    vector<Face>::size_type nf = faces.size();
    vector<vec>::size_type nv = vertices.size();

    for (vector<Face>::size_type i = 0; i < nf; i++) {
        const point &p0 = vertices[faces[i][0]];
        const point &p1 = vertices[faces[i][1]];
        const point &p2 = vertices[faces[i][2]];
        vec a = p0-p1, b = p1-p2, c = p2-p0;
        float l2a = len2(a), l2b = len2(b), l2c = len2(c);
        vec facenormal = a CROSS b;
        normals[faces[i][0]] += facenormal * (1.0f / (l2a * l2c));
        normals[faces[i][1]] += facenormal * (1.0f / (l2b * l2a));
        normals[faces[i][2]] += facenormal * (1.0f / (l2c * l2b));
    }

    for (vector<vec>::size_type i = 0; i < nv; i++)
        normalize(normals[i]);
}

// Compute per-vertex point areas
void Mesh::compute_pointareas()
{
	if (pointareas.size() == vertices.size())
		return;
    compute_faces();

    vector<Face>::size_type nf = faces.size();
    vector<point>::size_type nv = vertices.size();
	pointareas.clear();
	pointareas.resize(nv);
	cornerareas.clear();
	cornerareas.resize(nf);

	for (vector<Face>::size_type i = 0; i < nf; i++) {
		// Edges
		vec e[3] = { vertices[faces[i][2]] - vertices[faces[i][1]],
			     vertices[faces[i][0]] - vertices[faces[i][2]],
			     vertices[faces[i][1]] - vertices[faces[i][0]] };

		// Compute corner weights
		float area = 0.5f * len(e[0] CROSS e[1]);
		float l2[3] = { len2(e[0]), len2(e[1]), len2(e[2]) };
		float ew[3] = { l2[0] * (l2[1] + l2[2] - l2[0]),
				l2[1] * (l2[2] + l2[0] - l2[1]),
				l2[2] * (l2[0] + l2[1] - l2[2]) };
		if (ew[0] <= 0.0f) {
			cornerareas[i][1] = -0.25f * l2[2] * area /
					    (e[0] DOT e[2]);
			cornerareas[i][2] = -0.25f * l2[1] * area /
					    (e[0] DOT e[1]);
			cornerareas[i][0] = area - cornerareas[i][1] -
					    cornerareas[i][2];
		} else if (ew[1] <= 0.0f) {
			cornerareas[i][2] = -0.25f * l2[0] * area /
					    (e[1] DOT e[0]);
			cornerareas[i][0] = -0.25f * l2[2] * area /
					    (e[1] DOT e[2]);
			cornerareas[i][1] = area - cornerareas[i][2] -
					    cornerareas[i][0];
		} else if (ew[2] <= 0.0f) {
			cornerareas[i][0] = -0.25f * l2[1] * area /
					    (e[2] DOT e[1]);
			cornerareas[i][1] = -0.25f * l2[0] * area /
					    (e[2] DOT e[0]);
			cornerareas[i][2] = area - cornerareas[i][0] -
					    cornerareas[i][1];
		} else {
			float ewscale = 0.5f * area / (ew[0] + ew[1] + ew[2]);
			for (int j = 0; j < 3; j++)
				cornerareas[i][j] = ewscale * (ew[(j+1)%3] +
							       ew[(j+2)%3]);
		}
		pointareas[faces[i][0]] += cornerareas[i][0];
		pointareas[faces[i][1]] += cornerareas[i][1];
		pointareas[faces[i][2]] += cornerareas[i][2];
	}
}

void Mesh::compute_faces()
{
    // for next usage
    // to unpack triangle stripes
}

void Mesh::compute_adjacent_vertices()
{
    if (vertices.empty())
    {
        return;
    }
    if (adjacent.size() == vertices.size())
    {
        return;
    }

    for (unsigned int i=0; i<faces.size(); i++)
    {
        std::vector<int>::size_type v[3] = {faces[i][0], faces[i][1], faces[i][2]};

        adjacent[v[0]].insert(v[1]);
        adjacent[v[0]].insert(v[2]);

        adjacent[v[1]].insert(v[0]);
        adjacent[v[1]].insert(v[2]);

        adjacent[v[2]].insert(v[0]);
        adjacent[v[2]].insert(v[1]);
    }
}

__forceinline float Mesh::GaussianWeigtedAverage(unsigned int i, float delta)
{
    // all variables are the same with the equation in the paper.
    std::set<std::vector<int>::size_type> N;
    std::set<std::vector<int>::size_type> pN; // potential N;
    std::set<std::vector<int>::size_type> dN; // discarded N;

    pN.insert(i);

    // compute N
    while (!pN.empty())
    {
        std::vector<int>::size_type v = *(pN.begin());
        pN.erase(v);

        if (dist(vertices[i], vertices[v]) <= 2*delta)
        {
            N.insert(v);
        }
        else
        {
            dN.insert(v);
            continue; // if the vertex is out of delta range, then his neighbor is also discarded.
        }

        for (std::set<unsigned int>::iterator iter = adjacent[v].begin(); iter != adjacent[v].end(); ++iter)   
        {
            if (N.count(*iter) == 0 && dN.count(*iter) == 0)
            {
                pN.insert(*iter);
            }
        }
    }

    // compute the result
    double numerator=0;
    double denominator=0;
    for (std::set<unsigned int>::iterator iter = N.begin(); iter != N.end(); ++iter)
    {
        std::vector<unsigned int>::size_type idx = *iter;
        numerator += (curv1[idx] + curv2[idx])
            *exp(-dist2(vertices[i],vertices[idx])/(2*delta*delta));
        denominator += exp(-dist2(vertices[i],vertices[idx])/(2*delta*delta));
    }
    return float(numerator/denominator);
}

void Mesh::compute_mesh_saliency()
{
    if (vertices.empty())
    {
        return;
    }

    if (saliency.size() == vertices.size())
    {
        return;
    }

    // compute curvature;
    compute_curvatures();
    compute_bsphere();// bounding box
    compute_adjacent_vertices();

    epsilon = bsphere.r*2;
    // according to the paper , 0.3% of the diagonal of the BBox, here use diameter to approximate.
    epsilon *= 0.003f; 
    saliency.resize(vertices.size(),0);

    // 5 is the number of scales
    float delta[5] = {2*epsilon, 3*epsilon, 4*epsilon, 5*epsilon, 6*epsilon};
    // for each vertices, compute 
    for (unsigned int i=0; i<vertices.size(); i++)
    {
        float S[5];
        for (unsigned int j=0; j<5; j++)
        {
            S[j] = abs(GaussianWeigtedAverage(i,delta[j]) - GaussianWeigtedAverage(i,2*delta[j]));
        }
        // then, combine using a linear interpolation
        saliency[i] = (S[0]+S[1]+S[2]+S[3]+S[4])/5.0f;
    }

    // find max saliency and then, map to a color space
    // coz there always exist few bad curvature points, we just considering 99% of them
    float maxs = 0;
    for (unsigned int i=0; i<saliency.size(); i++)
    {
        if (saliency[i] > maxs)
        {
            maxs = saliency[i];
        }
    }
    int sz = int(maxs) + 2;

    unsigned int *histogram = new unsigned int[sz];
    memset(histogram,0,sz*sizeof(unsigned int));

    for (unsigned int i=0; i<saliency.size(); i++)
    {
        histogram[int(saliency[i]+0.5)] += 1;
    }

    // considering 99% vertices of the model
    const float PERC = 0.90f;
    int upperb=0;
    unsigned int count=histogram[0];
    while (float(count)/saliency.size() < PERC)
    {
        upperb++;
        count += histogram[upperb];
    }

    // map to color space
    // map saliency to color space
    // positive, red,
    // negative, green
    saliency_mapped_color.resize(saliency.size()*3);
    for (unsigned int i=0; i<saliency.size(); i++)
    {
        if (saliency[i] >= 0)
        {
            saliency_mapped_color[3*i + 0] = saliency[i]/upperb > 1 ? 1 : saliency[i]/upperb;
            saliency_mapped_color[3*i + 1] = 0;
            saliency_mapped_color[3*i + 2] = 0;
        }
        else
        {
            saliency_mapped_color[3*i + 0] = 0;
            saliency_mapped_color[3*i + 1] = (-saliency[i]/upperb > 1) ? 1 : (-saliency[i]/upperb);
            saliency_mapped_color[3*i + 2] = 0;
        }
    }

}

bool Mesh::read_segmentation(const wstring filename)
{
    if (filename.empty())
        return false;

    ifstream file;
    file.open(filename.c_str());
    if (file.fail()) {
        cout<<"file open error in read_segmentation()\n";
        return false;
    }

    unsigned int s;
    segmentation.clear();
    while(!file.eof())
    {
        file>>s;
        segmentation.push_back(s);
    }
    segmentation.pop_back(); //kick the eof last one

    file.close();
    return true;
}