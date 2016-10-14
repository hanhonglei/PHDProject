#include "com.h"

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <boost/random.hpp>

const float PI = 3.141592654f;


#define MIN(x, y)	( (x)<(y) ? (x) : (y) )
#define MAX(x, y)	( (x)>(y) ? (x) : (y) )
//
//HSB是以人类对颜色的感觉为依据而建立的，在HSB模式中，所有的颜色都是根据下列三种基本特征来进行形容和描述的。
//	色相（Hue）：是由物体发射或反射出来的颜色。它是根据色彩在一个0-360的标准色盘上的位置来决定的，通常以颜色的名称来辨识，
// 例如红、橙和绿等等，红色在0度，绿色在120度，蓝色在240度。它基本上是RGB模式全色度的饼状图。
//	饱和度（Saturation）：表示色彩的纯度，有时也被称为彩度，为0%时为灰色。白、黑和其他灰色色彩都没有饱和度的。在最大饱和度时，每一色相具有最纯的色光。
//	亮度（Brightness）：是指颜色相对的亮度和暗度，通常是用为0%（黑）到100%（白）的方式进行测定。
void Value2RGB(float v, float max, float min ,float rgb[3]) 
{  
	float R, G, B, H, S, I;
	H = (4/6.0f)*(max - v)/(max - min), S = 1, I = 1;

	if (I==0.0f) {
		// black image
		R = G = B = 0;
	}
	else {
		if (S==0.0f) {
			// grayscale image
			R = G = B = I;
		}
		else {
			float domainOffset = 0.0f;
			if (H<1.0f/6.0f) {	// red domain; green acending
				domainOffset = H;
				R = I;
				B = I * (1-S);
				G = B + (I-B)*domainOffset*6;
			}
			else {
				if (H<2.0f/6) {	// yellow domain; red acending
					domainOffset = H - 1.0f/6.0f;
					G = I;
					B = I * (1-S);
					R = G - (I-B)*domainOffset*6;
				}
				else {
					if (H<3.0f/6) {	// green domain; blue descending
						domainOffset = H-2.0f/6;
						G = I;
						R = I * (1-S);
						B = R + (I-R)*domainOffset * 6;
					}
					else {
						if (H<4.0f/6) {	// cyan domain, green acsending
							domainOffset = H - 3.0f/6;
							B = I;
							R = I * (1-S);
							G = B - (I-R) * domainOffset * 6;
						}
						else {
							if (H<5.0f/6) {	// blue domain, red ascending
								domainOffset = H - 4.0f/6;
								B = I;
								G = I * (1-S);
								R = G + (I-G) * domainOffset * 6;
							}
							else {	 // magenta domain, blue descending
								domainOffset = H - 5.0f/6;
								R = I;
								G = I * (1-S);
								B = R - (I-G) * domainOffset * 6;
							}
						}
					}
				}
			}
		}
	}

	rgb[0] = R;
	rgb[1] = G;
	rgb[2] = B;
}

void generate_viewpoint_candidates(size_t num, std::vector<point>& viewpoint_candidate_, int style)
{
    // clean and init
    viewpoint_candidate_.clear();
    viewpoint_candidate_.resize(num);

    time_t tm;
    srand((unsigned)time(&tm));  //set seed;

    if (style == -1)
    {
        // old style designed by myself, 
        for (size_t i=0; i<num; i++)
        {
            float alpha = (rand()%360)/360.0f*2*PI;        //angle with positive-x-axis
            float beta = (rand()%360)/360.0f*2*PI - PI;  //angle with xy-plane

            viewpoint_candidate_[i][0] = cos(beta)*sin(alpha);
            viewpoint_candidate_[i][1] = cos(beta)*cos(alpha);
            viewpoint_candidate_[i][2] = sin(beta);
        }
        //for (unsigned int i=0; i<num; i++)
        //{
        //    float alpha = i*30.0/180.0*PI; 
        //    float beta = 0.0;

        //    viewpoint_candidate_[i][0] = cos(beta)*sin(alpha);
        //    viewpoint_candidate_[i][1] = sin(beta);
        //    viewpoint_candidate_[i][2] = cos(beta)*cos(alpha);
        //}
    }
    else if(style == 0)
    {
        // method in harmonic lighting tutorial
        // generation of two series of random number, which are independent
        // when sampling ,use it to multiple bounding sphere radius
        boost::mt19937 rng_a;                 // produces randomness out of thin air
        boost::uniform_01<> dist_a;
        boost::variate_generator<boost::mt19937&, boost::uniform_01<> >
            rand_a(rng_a, dist_a);             // glues randomness with mapping
        rng_a.seed(rand());
        //rng_a.seed(1);
        std::vector<float> ra;
        for (unsigned int i=0; i<num; i++)
        {
            ra.push_back((float)rand_a());
        }

        boost::mt11213b rng_b;                 // produces randomness out of thin air
        boost::uniform_01<> dist_b;
        boost::variate_generator<boost::mt11213b&, boost::uniform_01<> >
            rand_b(rng_b, dist_b);             // glues randomness with mapping
        rng_b.seed(rand());
        //rng_b.seed(8);
        std::vector<float> rb;
        for (unsigned int i=0; i<num; i++)
        {
            rb.push_back((float)rand_b());
        }
        float theta = 0;
        float phi = 0;
        for (unsigned int i=0; i<num; i++)
        {
            theta = 2*acos(sqrt(1-ra[i]));
            phi = 2*PI*rb[i];

            viewpoint_candidate_[i][0] = sin(theta)*cos(phi);
            viewpoint_candidate_[i][1] = sin(theta)*sin(phi);
            viewpoint_candidate_[i][2] = cos(theta);
        }
    }
}


