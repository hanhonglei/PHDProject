#include "ppm.h"

bool PPM::load(std::wstring fn)
{
    std::ifstream f(fn.c_str(), std::ios::binary|std::ios::in);
    f>>version;
    if(version != std::string("P6")){
        std::cout<<"Unsupported PPM file format... Only support P6"<<std::endl;
        return false;
    }
    f>>width>>height;
    int ncolor;
    f>>ncolor;
    if (!(ncolor>0 && ncolor < 256))
    {
        std::cout<<"Unsupported PPM file format... Only support 256 colors"<<std::endl;
        return false;
    }
    f.seekg(1,std::ios::cur);   // the Return after 255, following is the color data

    // if there is already data in the object.
    if (data!=0)
    {
        delete [] data;
    }
    data = new unsigned char[width*height*3];
    f.read((char*)data,width*height*3);

    // swap i and height-i rows, i=0, 1, ... [height/2]
	// top-to-bottom
	for (unsigned int i = 0; i < height/2; i++) {
		unsigned char *row1 = data + 3 * width * i;
		unsigned char *row2 = data + 3 * width * (height - 1 - i);
		for (unsigned int j = 0; j < 3 * width; j++)
            std::swap(row1[j], row2[j]);
	}

    return true;
}

bool PPM::save(std::wstring fn)
{
    std::ofstream f(fn.c_str(), std::ios::binary|std::ios::out);
    f<<version<<std::endl;
    f<<width<<' '<<height<<std::endl;
    f<<255<<std::endl;

    // swap i and height-i rows, i=0, 1, ... [height/2]
	// top-to-bottom
	for (unsigned int i = 0; i < height/2; i++) {
		unsigned char *row1 = data + 3 * width * i;
		unsigned char *row2 = data + 3 * width * (height - 1 - i);
		for (unsigned int j = 0; j < 3 * width; j++)
            std::swap(row1[j], row2[j]);
	}

    f.write((char*)data,width*height*3*sizeof(unsigned char));

    // swap i and height-i rows, i=0, 1, ... [height/2]
	// top-to-bottom
	for (unsigned int i = 0; i < height/2; i++) {
		unsigned char *row1 = data + 3 * width * i;
		unsigned char *row2 = data + 3 * width * (height - 1 - i);
		for (unsigned int j = 0; j < 3 * width; j++)
            std::swap(row1[j], row2[j]);
	}
	f.close();
    return true;
}