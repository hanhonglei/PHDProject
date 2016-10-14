#pragma  once

#include <iostream>
#include <string>
#include <fstream>
#include <boost/smart_ptr.hpp>
#include "mesh.h"

// 3d volume data
class VolumeData{
public:
    static const unsigned int GRANULARITY = 256*3; // the voxelization granularity, approximating x + y + z

    // size along x, y, z axis
    int sx;
    int sy;
    int sz;

    unsigned int* data;

    VolumeData(){data=0;}
    ~VolumeData(){if(data != 0) delete [] data;}

    bool load(std::wstring fn);
    bool voxelize(Mesh *mesh);
    bool save(std::wstring fn);
};