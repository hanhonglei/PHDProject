#ifndef MXSTDSLIM_INCLUDED // -*- C++ -*-
#define MXSTDSLIM_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Core simplification interface.  The MxStdSlim class defines the
  interface which all simplification classes conform to.

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxStdSlim.h,v 1.4 1998/11/19 01:57:34 garland Exp $

 ************************************************************************/

#include "MxStdModel.h"
#include "MxHeap.h"

#define MX_PLACE_ENDPOINTS 0
#define MX_PLACE_ENDORMID  1
#define MX_PLACE_LINE      2
#define MX_PLACE_OPTIMAL   3

#define MX_WEIGHT_UNIFORM       0
#define MX_WEIGHT_AREA          1
#define MX_WEIGHT_ANGLE         2
// 使用平均Q，用来决定最大几何误差。如果使用原始的递加Q的话，得到的误差是当前顶点离相关平面的距离总和 [5/21/2012 Han]
// 每个Q的结构体中保存一个int型成员，保存Q所代表的多边形个数
#define MX_WEIGHT_AVERAGE       3
#define MX_WEIGHT_AREA_AVG      4
#define MX_WEIGHT_RAWNORMALS    5

class MxStdSlim
{
protected:
    MxStdModel *m;
    MxHeap heap;

public:
    unsigned int valid_verts;
    unsigned int valid_faces;
    bool is_initialized;

    int placement_policy;
    int weighting_policy;
    bool will_join_only;

    double boundary_weight;
    double compactness_ratio;
    double meshing_penalty;
    double local_validity_threshold;
    uint vertex_degree_limit;

public:
    MxStdSlim(MxStdModel *m0);

    virtual void initialize() = 0;
    virtual bool decimate(uint) = 0;

    MxStdModel& model() { return *m; }
	// 利用距离作为衡量因素进行边折叠计算 [3/25/2012 Han]
	// 计算依据:计算这个距离下,一个像素反投影到模型所在位置的距离平方QError。进行边折叠简化，直到误差超过QError
	virtual double AdjustLOD(float d, float r, double l) = 0;
	// QError is the quadric distance of collapsed vertex to it's original position
	virtual bool DecimateByError(double QError) = 0;


};

// MXSTDSLIM_INCLUDED
#endif
