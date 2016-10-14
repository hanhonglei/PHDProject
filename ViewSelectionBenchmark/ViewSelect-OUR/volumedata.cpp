#include "volumedata.h"
//#include "debug.h"


bool VolumeData::load(std::wstring fn)
{
    Mesh *mesh = Mesh::read(fn);
    //ASSERT(mesh);

    voxelize(mesh);

    delete mesh;
    return true;
}

bool VolumeData::voxelize(Mesh *mesh)
{
    //ASSERT(mesh);
    mesh->compute_bsphere();

    float lx = mesh->bsphere.center[0] - mesh->bsphere.r; // low bounder of x
    float ux = mesh->bsphere.center[0] + mesh->bsphere.r; // up bounder of x
    float dx = ux - lx; // length of the bounding box along x axis
    float ly = mesh->bsphere.center[1] - mesh->bsphere.r;
    float uy = mesh->bsphere.center[1] + mesh->bsphere.r;
    float dy = uy - ly;
    float lz = mesh->bsphere.center[2] - mesh->bsphere.r;
    float uz = mesh->bsphere.center[2] + mesh->bsphere.r;
    float dz = uz - lz;

    float vs = (dx+dy+dz)/GRANULARITY;  // voxel size
    sx = unsigned int(dx / vs);         // size of the whole data along x axis
    sy = unsigned int(dy / vs);
    sz = unsigned int(dz / vs);
    data = new unsigned int[sz * sy * sz];
    //ASSERT(data);
    memset(data,0,sz*sy*sz*sizeof(unsigned int));

    for (unsigned int i=0; i<mesh->vertices.size(); i++)
    {
        unsigned int ix = unsigned int((mesh->vertices[i][0] - lx) / vs);
        unsigned int iy = unsigned int((mesh->vertices[i][1] - ly) / vs);
        unsigned int iz = unsigned int((mesh->vertices[i][2] - lz) / vs);
        
        data[ix + iy * sx + iz * sx * sy] += 1;
    }
    return true;
}

bool VolumeData::save(std::wstring fn)
{
    return true;
}