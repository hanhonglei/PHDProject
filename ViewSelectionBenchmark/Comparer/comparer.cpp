#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "D:\Projects\ViewSelectionBenchmark\HumanTest\Vec.h"


const unsigned int TM = 45;

void compare2human(const std::string & src1, const std::string & src2, const std::string & ofname)
{
    std::ifstream f1(src1.c_str());
    std::ifstream f2(src2.c_str());
    std::ofstream of(ofname.c_str());
    std::ofstream oftime((ofname+std::string(".time")).c_str());

    if (!f1.good() || !f2.good() || !of.good())
    {
        std::cout<<"File open error in: "<<src1<<" "<<src2<<" "<<ofname<<std::endl;
        return;
    }

    float tc; //time cost
    point v1;
    point v2;
    std::string linebuffer;
    std::stringstream streambuffer;
    for (unsigned int i=0; i<45; ++i)
    {
        // parse human viewpoint position of Model-i;
        getline(f1,linebuffer);
        streambuffer.clear();
        streambuffer.str(linebuffer);
        streambuffer>>v1;
        getline(f2,linebuffer);
        streambuffer.clear();
        streambuffer.str(linebuffer);
        streambuffer>>tc>>v2;

        // compute the angle 
        normalize(v1);
        normalize(v2);
        float ang = v1 DOT v2;
        ang = ang>1?1:ang;
        ang = ang<0?0:ang;
        ang = acos(ang);

        oftime<<tc<<std::endl;
        of<<ang<<std::endl;
    }
    f1.close();
    f2.close();
    of.close();
}

int main(int argc, char **argv)
{
    std::string humanresult = "humanresult";

    std::string ourresult = "ourresult";
    std::string veresult = "veresult";
    std::string vmiresult = "vmiresult";
    std::string msresult = "msresult";
    std::string dpresult = "dpresult";

    std::string human_our = "human-our";
    compare2human(humanresult,ourresult,human_our);

    std::string human_ve = "human-ve";
    compare2human(humanresult,veresult,human_ve);

    std::string human_vmi = "human-vmi";
    compare2human(humanresult,vmiresult,human_vmi);

    std::string human_ms = "human-ms";
    compare2human(humanresult,msresult,human_ms);

    std::string human_dp = "human-dp";
    compare2human(humanresult,dpresult,human_dp);
    return 0;
}