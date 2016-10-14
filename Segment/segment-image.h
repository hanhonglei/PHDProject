/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include <cstdlib>
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "segment-graph.h"
#include "filter.h"
#include <stdmix.h>
#include <MxTimer.h>
#include "qslim.h"
#include <list>

// random color
rgb random_rgb();

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
			 int x1, int y1, int x2, int y2) ;

void  segment_image(MxStdModel *pMesh, float sigma, float c, int min_size,int *num_ccs, bool bFace);
void  segment_image_face(MxStdModel *pMesh, float sigma, float c, int min_size,int *num_ccs);
void  segment_image_vert(MxStdModel *pMesh, float sigma, float c, int min_size,int *num_ccs);


#endif