#pragma once

#include <iostream>
#include <string>
#include <fstream>

class PPM{
public:
    std::string version;
    unsigned int width;
    unsigned int height;
    unsigned char *data;

    PPM(){data=0;}
    ~PPM(){if(data != 0) delete [] data;}

    bool load(std::wstring fn);
    bool save(std::wstring fn);
};
