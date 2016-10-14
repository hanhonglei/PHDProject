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

#include <cstdio>
#include <cstdlib>
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "segment-image.h"

//////////////////////////////////////////////////////////////////////////
// how to use segment for image
/*
int main(int argc, char **argv) {
if (argc != 6) {
fprintf(stderr, "usage: %s sigma k min input(ppm) output(ppm)\n", argv[0]);
getchar();
return 1;
}

float sigma = atof(argv[1]);
float k = atof(argv[2]);
int min_size = atoi(argv[3]);

printf("loading input image.\n");
image<rgb> *input = loadPPM(argv[4]);

printf("processing\n");
int num_ccs; 
image<rgb> *seg = segment_image(input, sigma, k, min_size, &num_ccs); 
savePPM(seg, argv[5]);

printf("got %d components\n", num_ccs);
printf("done! uff...thats hard work.\n");
getchar();

return 0;
}

*/


// random color
rgb random_rgb(){ 
	rgb c;
	double r;

	c.r = (uchar)rand();
	c.g = (uchar)rand();
	c.b = (uchar)rand();

	return c;
}

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
	int x1, int y1, int x2, int y2) {
		return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
			square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
			square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}

/*
* Segment an image
*
* Returns a color image representing the segmentation.
*
* im: image to segment.
* sigma: to smooth the image.
* c: constant for treshold function.
* min_size: minimum component size (enforced by post-processing stage).
* num_ccs: number of connected components in the segmentation.
// */
//image<rgb> *segment_image(image<rgb> *im, float sigma, float c, int min_size,
//			  int *num_ccs) {
//  int width = im->width();
//  int height = im->height();
//
//  image<float> *r = new image<float>(width, height);
//  image<float> *g = new image<float>(width, height);
//  image<float> *b = new image<float>(width, height);
//
//  // smooth each color channel  
//  for (int y = 0; y < height; y++) {
//    for (int x = 0; x < width; x++) {
//      imRef(r, x, y) = imRef(im, x, y).r;
//      imRef(g, x, y) = imRef(im, x, y).g;
//      imRef(b, x, y) = imRef(im, x, y).b;
//    }
//  }
//  image<float> *smooth_r = smooth(r, sigma);
//  image<float> *smooth_g = smooth(g, sigma);
//  image<float> *smooth_b = smooth(b, sigma);
//  delete r;
//  delete g;
//  delete b;
// 
//  // build graph
//  edge *edges = new edge[width*height*4];
//  int num = 0;
//  for (int y = 0; y < height; y++) {
//    for (int x = 0; x < width; x++) {
//      if (x < width-1) {
//	edges[num].a = y * width + x;
//	edges[num].b = y * width + (x+1);
//	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
//	num++;
//      }
//
//      if (y < height-1) {
//	edges[num].a = y * width + x;
//	edges[num].b = (y+1) * width + x;
//	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
//	num++;
//      }
//
//      if ((x < width-1) && (y < height-1)) {
//	edges[num].a = y * width + x;
//	edges[num].b = (y+1) * width + (x+1);
//	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
//	num++;
//      }
//
//      if ((x < width-1) && (y > 0)) {
//	edges[num].a = y * width + x;
//	edges[num].b = (y-1) * width + (x+1);
//	edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
//	num++;
//      }
//    }
//  }
//  delete smooth_r;
//  delete smooth_g;
//  delete smooth_b;
//
//  // segment
//  universe *u = segment_graph(width*height, num, edges, c);
//  
//  // post process small components
//  for (int i = 0; i < num; i++) {
//    int a = u->find(edges[i].a);
//    int b = u->find(edges[i].b);
//    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
//      u->join(a, b);
//  }
//  delete [] edges;
//  *num_ccs = u->num_sets();
//
//  image<rgb> *output = new image<rgb>(width, height);
//
//  // pick random colors for each component
//  rgb *colors = new rgb[width*height];
//  for (int i = 0; i < width*height; i++)
//    colors[i] = random_rgb();
//  
//  for (int y = 0; y < height; y++) {
//    for (int x = 0; x < width; x++) {
//      int comp = u->find(y * width + x);
//      imRef(output, x, y) = colors[comp];
//    }
//  }  
//
//  delete [] colors;  
//  delete u;
//
//  return output;
//}
/*
* Segment an image
*
* Returns a color image representing the segmentation.
*
* im: image to segment.
* sigma: to smooth the image.
* c: constant for treshold function.
* min_size: minimum component size (enforced by post-processing stage).
* num_ccs: number of connected components in the segmentation.
// */

// 返回值代表分块的数目 [7/10/2012 Han]
void  segment_image(MxStdModel *pMesh, float sigma, float c, int min_size,int *num_ccs, bool bFace)
{
	if (bFace)
		segment_image_face(pMesh, sigma, c, min_size,num_ccs);
	else
		segment_image_vert(pMesh, sigma, c, min_size,num_ccs);
}

// 使用面片作为顶点、相邻面片看成边进行聚类 [7/10/2012 Han]
void  segment_image_face(MxStdModel *pMesh, float sigma, float c, int min_size,int *num_ccs)
{
	std::list<edge> edges;
	std::map<unsigned int, unsigned int> oldID_newID;

	int vertNum = 0;
	for (MxFaceID f = 0; f < pMesh->face_count(); f++)
	{
		oldID_newID.insert(make_pair(f, vertNum));
		if (!pMesh->face_is_valid(f))
			continue;

		float cfi = 0.0f;
		//for(MxVertexID v = 0; v < 3; v++)
		//	cfi += pMesh->vertex(pMesh->face(f)[v]).view_importance;
		cfi = pMesh->face(f).view_importance;

		std::map<int ,float> neighbours;

		// 每个顶点的相邻三角形都算作本三角形的邻居 [7/11/2012 Han]
		//for(MxVertexID v = 0; v < 3; v++)
		//{
		//	const MxFaceList& nFL = pMesh->neighbors(pMesh->face(f)[v]);
		//	for(int fn = 0; fn < nFL.length(); fn++)
		//	{
		//		if (!pMesh->face_is_valid(nFL(fn)))
		//			continue;

		//		neighbours.insert(make_pair(nFL[fn], pMesh->face(nFL(fn)).view_importance));
		//	}
		//}
		// 只有边相邻才算作邻居 [7/11/2012 Han]
		for(MxVertexID v = 0; v < 3; v++)
		{
			MxFaceList edgeN;
			edgeN.reset();
			pMesh->collect_edge_neighbors(pMesh->face(f)[v], pMesh->face(f)[(v+1)%3], edgeN);
			for(int fn = 0; fn < edgeN.length(); fn++)
			{
				neighbours.insert(make_pair(edgeN[fn], pMesh->face(edgeN(fn)).view_importance));
				if (edgeN(fn) >= pMesh->face_count() || !pMesh->face_is_valid(edgeN(fn)))
				{
					int test = 0;
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////

		std::map<int, float>::iterator  It = neighbours.begin();
		while(It != neighbours.end())
		{
			//const MxFace &fS = pMesh->face(nFL[fn]);

			if (It->first > f && pMesh->face_is_valid(It->first))
			{
				edge eg;
				eg.a = f;
				eg.b = It->first;
				float nfi = 0.0f;
				//for(MxVertexID v = 0; v < 3; v++)
				//	nfi += pMesh->vertex(fS[v]).view_importance;
				nfi = It->second;
				eg.w = abs(cfi - nfi);
				edges.push_back(eg);
			}
			It++;
		}
		vertNum++;
	}

	// segment
	universe *u = segment_graph(/*pMesh->face_count()*/vertNum, edges.size(), edges, c, oldID_newID);		// 使用原始的三角形数目，否则报错

	// post process small components
	std::list<edge>::iterator it;

	for (it=edges.begin(); it!=edges.end(); ++it){
		int newa = oldID_newID[it->a];
		int newb = oldID_newID[it->b];

		int a = u->find(newa/*it->a*/);
		int b = u->find(newb/*it->b*/);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}
	//delete [] edges;
	edges.clear();
	pMesh->nNumSeg = *num_ccs = u->num_sets();

	int nmax = 0, nmin = INT_MAX;


	for (MxFaceID f = 0; f < pMesh->face_count(); f++)
	{
		if (!pMesh->face_is_valid(f))
			continue;
		unsigned int newf = oldID_newID[f];
		int comp = u->find(newf);

		pMesh->face(f).from = comp;

		for(MxVertexID v = 0; v < 3; v++)
			pMesh->vertex(pMesh->face(f)[v]).from = comp;


		if (comp > nmax)
			nmax = comp;
		if (comp < nmin)
			nmin = comp;			
	}

	delete u;
}

void  segment_image_vert(MxStdModel *pMesh, float sigma, float c, int min_size,int *num_ccs)
{
	std::list<edge> edges;
	std::map <unsigned int, unsigned int> oldID_newID;
	int vertNum = 0;
	for (MxVertexID i = 0; i < pMesh->vert_count(); i++)
	{
		oldID_newID.insert(make_pair(i, vertNum));
		if (!pMesh->vertex_is_valid(i))
			continue;

		MxVertexList star;
		star.reset();
		if (i == 13)
		{
			int test = 1;
			//  [2/17/2014 Han] test
			//const MxFaceList& N = 		pMesh->neighbors(i);
			//for(uint tttt=0; tttt<N.length(); tttt++)
			//{
			//	char buffer[64];
			//	sprintf(buffer, "F:\t%d\n", N(tttt));
			//	OutputDebugString(buffer);
			//	for(uint j=0; j<3; j++)
			//	{
			//		char buffer[64];
			//		int test = 1;
			//		sprintf(buffer, "V:%d\n",  pMesh->face(N(tttt))(j));
			//		OutputDebugString(buffer);

			//	}
			//}
			// end test [2/17/2014 Han]

		}
		
		pMesh->collect_vertex_star(i, star);
		//MxEdgeList h;
		//edge_around_vertex_circular(v, h);
		//for(uint j=0; j<h.length(); j++)
		for(int j = 0; j < star.length(); j++)
		{
			if (star(j) > i)
			{
				edge eg;
				eg.a = i;
				eg.b = star(j);
				eg.w = abs(pMesh->vertex(i).view_importance-pMesh->vertex(star(j)).view_importance);
				edges.push_back(eg);
			}
		}
		vertNum++;
				//char buffer[64];
				//sprintf(buffer, "%d\n", i);
				//OutputDebugString(buffer);
	}
	//for (int y = 0; y < height; y++) {
	//	for (int x = 0; x < width; x++) {
	//		if (x < width-1) {
	//			edges[num].a = y * width + x;
	//			edges[num].b = y * width + (x+1);
	//			edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
	//			num++;
	//		}

	//		if (y < height-1) {
	//			edges[num].a = y * width + x;
	//			edges[num].b = (y+1) * width + x;
	//			edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
	//			num++;
	//		}

	//		if ((x < width-1) && (y < height-1)) {
	//			edges[num].a = y * width + x;
	//			edges[num].b = (y+1) * width + (x+1);
	//			edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
	//			num++;
	//		}

	//		if ((x < width-1) && (y > 0)) {
	//			edges[num].a = y * width + x;
	//			edges[num].b = (y-1) * width + (x+1);
	//			edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
	//			num++;
	//		}
	//	}
	//}
	//delete smooth_r;
	//delete smooth_g;
	//delete smooth_b;

	// segment
	universe *u = segment_graph(pMesh->vert_count()/*vertNum*/, edges.size(),  edges, c, oldID_newID);	// 使用原来的顶点数，否则空间不够，报错

	// post process small components
	std::list<edge>::iterator it;

	for (it=edges.begin(); it!=edges.end(); ++it){
		int newa = oldID_newID[it->a];
		int newb = oldID_newID[it->b];

		int a = u->find(it->a);
		int b = u->find(it->b);
		if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
			u->join(a, b);
	}
	//delete [] edges;
	edges.clear();
	pMesh->nNumSeg = *num_ccs = u->num_sets();

	//image<rgb> *output = new image<rgb>(width, height);

	// pick random colors for each component
	//rgb *colors = new rgb[width*height];
	//for (int i = 0; i < width*height; i++)
	//	colors[i] = random_rgb();

	//for (int y = 0; y < height; y++) {
	//	for (int x = 0; x < width; x++) {

	int nmax = 0, nmin = INT_MAX;
	//int nMaxNum = pMesh->vert_count();
	//int *compNum = new int[nMaxNum];
	//memset(compNum, 0, sizeof(int)*nMaxNum);
	//float *compSaliency = new float[nMaxNum];
	//memset(compSaliency, 0.f, sizeof(float)*nMaxNum);
	for (MxVertexID i = 0; i < pMesh->vert_count(); i++){
		if (!pMesh->vertex_is_valid(i))
			continue;
		int newi = oldID_newID[i];
		int comp = u->find(newi);
		//imRef(output, x, y) = colors[comp];
		// test 将得到的类别设置为importance [7/10/2012 Han]
		// 借用from指定当前顶点所属的分块 [7/10/2012 Han]
		pMesh->vertex(i).from = comp;
		//compSaliency[comp] += pMesh->vertex(i).view_importance ;
		//compNum[comp]++;
		if (comp > nmax)
			nmax = comp;
		if (comp < nmin)
			nmin = comp;			
	}
	// 依据顶点的归属，计算面片的归属 [7/12/2012 Han]
	int compsVerts[3];
	for (MxFaceID f = 0; f < pMesh->face_count(); f++)
	{
		if (!pMesh->face_is_valid(f))
			continue;

		for(MxVertexID v = 0; v < 3; v++)
			compsVerts[v] = pMesh->vertex(pMesh->face(f)[v]).from ;
		int comp = compsVerts[0];
		// 选最多的 [7/12/2012 Han]
		comp = (compsVerts[1] == compsVerts[2])?compsVerts[1]:( (compsVerts[0] == compsVerts[1])?compsVerts[0]:compsVerts[2]);

		pMesh->face(f).from = comp;

		if (comp > nmax)
			nmax = comp;
		if (comp < nmin)
			nmin = comp;			
	}

	// 依照分组，对每个顶点进行重新复制 [7/10/2012 Han]
	// 暂时不做，放到外部，按照需要再做
	//for (int i = 0; i < nMaxNum; i++)
	//	if (compNum[i] != 0)
	//		compSaliency[i] /= float(compNum[i]);

	//for (MxVertexID i = 0; i < pMesh->vert_count(); i++)
	//	pMesh->vertex(i).view_importance = compSaliency[pMesh->vertex(i).from];


	//delete [] colors;  
	delete u;

	//return output;
}

bool operator<(const edge &a, const edge &b) {
	return a.w < b.w;
}

bool compar(const edge &a, const edge &b) {
	return a.w < b.w;
}
/*
* Segment a graph
*
* Returns a disjoint-set forest representing the segmentation.
*
* num_vertices: number of vertices in graph.
* num_edges: number of edges in graph
* edges: array of edges.
* c: constant for treshold function.
*/
//universe *segment_graph(int num_vertices, int num_edges, edge *edges, 
//			float c) { 
universe *segment_graph(int num_vertices, int num_edges, std::list<edge> &edges, float c, std::map<unsigned int, unsigned int> &olo_newID) { 

	// sort edges by weight
	//std::sort(edges, edges + num_edges);
	edges.sort(compar);

	// make a disjoint-set forest
	universe *u = new universe(num_vertices);

	// init thresholds
	float *threshold = new float[num_vertices];
	for (int i = 0; i < num_vertices; i++)
		threshold[i] = THRESHOLD(1,c);

	// for each edge, in non-decreasing weight order...
	//for (int i = 0; i < num_edges; i++) {
	std::list<edge>::iterator it;

	for (it=edges.begin(); it!=edges.end(); ++it){


		//edge *pedge = &edges[i];
		int newa = olo_newID[it->a];
		int newb = olo_newID[it->b];

		// components conected by this edge
		int a = u->find(newa/*it->a*/);
		int b = u->find(newb/*it->b*/);
		if (a != b) {
			if ((it->w <= threshold[a]) &&
				(it->w <= threshold[b])) {
					u->join(a, b);
					a = u->find(a);
					threshold[a] = it->w + THRESHOLD(u->size(a), c);
			}
		}
	}

	// free up
	delete []threshold;
	return u;
}


static void normalize(std::vector<float> &mask) {
	int len = mask.size();
	float sum = 0;
	for (int i = 1; i < len; i++) {
		sum += fabs(mask[i]);
	}
	sum = 2*sum + fabs(mask[0]);
	for (int i = 0; i < len; i++) {
		mask[i] /= sum;
	}
}

/* make filters */
//#define MAKE_FILTER(name, fun)                                \
//	static std::vector<float> make_ ## name (float sigma) {       \
//	sigma = std::max(sigma, 0.01F);			      \
//	int len = (int)ceil(sigma * WIDTH) + 1;                     \
//	std::vector<float> mask(len);                               \
//	for (int i = 0; i < len; i++) {                             \
//	mask[i] = fun;                                            \
//	}                                                           \
//	return mask;                                                \
//}

//MAKE_FILTER(fgauss, exp(-0.5*square(i/sigma)));

/* convolve image with gaussian filter */
static image<float> *smooth(image<float> *src, float sigma) {
	std::vector<float> mask = make_fgauss(sigma);
	normalize(mask);

	image<float> *tmp = new image<float>(src->height(), src->width(), false);
	image<float> *dst = new image<float>(src->width(), src->height(), false);
	convolve_even(src, tmp, mask);
	convolve_even(tmp, dst, mask);

	delete tmp;
	return dst;
}

/* convolve image with gaussian filter */
image<float> *smooth(image<uchar> *src, float sigma) {
	image<float> *tmp = imageUCHARtoFLOAT(src);
	image<float> *dst = smooth(tmp, sigma);
	delete tmp;
	return dst;
}

/* compute laplacian */
static image<float> *laplacian(image<float> *src) {
	int width = src->width();
	int height = src->height();
	image<float> *dst = new image<float>(width, height);  

	for (int y = 1; y < height-1; y++) {
		for (int x = 1; x < width-1; x++) {
			float d2x = imRef(src, x-1, y) + imRef(src, x+1, y) -
				2*imRef(src, x, y);
			float d2y = imRef(src, x, y-1) + imRef(src, x, y+1) -
				2*imRef(src, x, y);
			imRef(dst, x, y) = d2x + d2y;
		}
	}
	return dst;
}


universe::universe(int elements) {
	elts = new uni_elt[elements];
	num = elements;
	for (int i = 0; i < elements; i++) {
		elts[i].rank = 0;
		elts[i].size = 1;
		elts[i].p = i;
	}
}

universe::~universe() {
	delete [] elts;
}

int universe::find(int x) {
	int y = x;
	while (y != elts[y].p)
		y = elts[y].p;
	elts[x].p = y;
	return y;
}

void universe::join(int x, int y) {
	if (elts[x].rank > elts[y].rank) {
		elts[y].p = x;
		elts[x].size += elts[y].size;
	} else {
		elts[x].p = y;
		elts[y].size += elts[x].size;
		if (elts[x].rank == elts[y].rank)
			elts[y].rank++;
	}
	num--;
}