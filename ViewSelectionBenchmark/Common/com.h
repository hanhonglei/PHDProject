# pragma once

#include <vector>
#include "Vec.h"

void generate_viewpoint_candidates(size_t num, std::vector<point>& viewpoint_candidate_, int style=0);
void Value2RGB(float v, float max, float min ,float rgb[3]) ;

