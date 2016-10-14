#include "glmatrix.h"
#include <cmath>
#include <cassert>
#define myPi acos(-1.0f)

//All this matrix are in row-major mode, for Cg usage

void buildIdentityMatrix(float m[16])
{
    m[0]  = 1;  m[1]  = 0;  m[2]  = 0;  m[3]  = 0;
    m[4]  = 0;  m[5]  = 1;  m[6]  = 0;  m[7]  = 0;
    m[8]  = 0;  m[9]  = 0;  m[10] = 1;  m[11] = 0;
    m[12] = 0;  m[13] = 0;  m[14] = 0;  m[15] = 1;
}

void buildOrthoMatrix(double left, double right,
                      double bottom, double top,
                      double near, double far,
                      float m[16])
{
#pragma warning(push)
#pragma warning(disable : 4244)
    m[0*4+0] = 2.0/(right-left);
    m[0*4+1] = 0.0;
    m[0*4+2] = 0.0;
    m[0*4+3] = -(right+left)/(right-left);

    m[1*4+0] = 0.0;
    m[1*4+1] = 2.0/(top-bottom);
    m[1*4+2] = 0.0;
    m[1*4+3] = -(top+bottom)/(top-bottom);

    m[2*4+0] = 0.0;
    m[2*4+1] = 0.0;
    m[2*4+2] = -2.0/(far-near);
    m[2*4+3] = -(far+near)/(far-near);

    m[3*4+0] = 0.0;
    m[3*4+1] = 0.0;
    m[3*4+2] = 0.0;
    m[3*4+3] = 1.0;
#pragma warning(pop)
}

void buildPerspectiveMatrix(double fieldOfView,
                            double aspectRatio,
                            double zNear, double zFar,
                            float m[16])
{
    double sine, cotangent, deltaZ;
    double radians = fieldOfView / 2.0 * myPi / 180.0;

    deltaZ = zFar - zNear;
    sine = sin(radians);
    // Should be non-zero to avoid division by zero.
    assert(deltaZ);
    assert(sine);
    assert(aspectRatio);
    cotangent = cos(radians) / sine;

#pragma warning(push)
#pragma warning(disable : 4244)
    m[0*4+0] = cotangent / aspectRatio;
    m[0*4+1] = 0.0;
    m[0*4+2] = 0.0;
    m[0*4+3] = 0.0;

    m[1*4+0] = 0.0;
    m[1*4+1] = cotangent;
    m[1*4+2] = 0.0;
    m[1*4+3] = 0.0;

    m[2*4+0] = 0.0;
    m[2*4+1] = 0.0;
    m[2*4+2] = -(zFar + zNear) / deltaZ;
    m[2*4+3] = -2 * zNear * zFar / deltaZ;

    m[3*4+0] = 0.0;
    m[3*4+1] = 0.0;
    m[3*4+2] = -1;
    m[3*4+3] = 0;
#pragma warning(pop)
}

void buildLookAtMatrix(double eyex, double eyey, double eyez,
                       double centerx, double centery, double centerz,
                       double upx, double upy, double upz,
                       float m[16])
{
    double x[3], y[3], z[3], mag;

    z[0] = eyex - centerx;
    z[1] = eyey - centery;
    z[2] = eyez - centerz;

    mag = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
    if (mag) {
        z[0] /= mag;
        z[1] /= mag;
        z[2] /= mag;
    }

    y[0] = upx;
    y[1] = upy;
    y[2] = upz;

    x[0] =  y[1]*z[2] - y[2]*z[1];
    x[1] = -y[0]*z[2] + y[2]*z[0];
    x[2] =  y[0]*z[1] - y[1]*z[0];

    y[0] =  z[1]*x[2] - z[2]*x[1];
    y[1] = -z[0]*x[2] + z[2]*x[0];
    y[2] =  z[0]*x[1] - z[1]*x[0];

    mag = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    if (mag) {
        x[0] /= mag;
        x[1] /= mag;
        x[2] /= mag;
    }

    mag = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
    if (mag) {
        y[0] /= mag;
        y[1] /= mag;
        y[2] /= mag;
    }

#pragma warning(push)
#pragma warning(disable : 4244)
    m[0*4+0] = x[0];  m[0*4+1] = x[1];
    m[0*4+2] = x[2];  m[0*4+3] = -x[0]*eyex + -x[1]*eyey + -x[2]*eyez;

    m[1*4+0] = y[0];  m[1*4+1] = y[1];
    m[1*4+2] = y[2];  m[1*4+3] = -y[0]*eyex + -y[1]*eyey + -y[2]*eyez;

    m[2*4+0] = z[0];  m[2*4+1] = z[1];
    m[2*4+2] = z[2];  m[2*4+3] = -z[0]*eyex + -z[1]*eyey + -z[2]*eyez;

    m[3*4+0] = 0.0;   m[3*4+1] = 0.0;  m[3*4+2] = 0.0;  m[3*4+3] = 1.0;
#pragma warning(pop)
}

void buildRotateMatrix(float angle,
                      float ax, float ay, float az,
                      float m[16])
{
    float radians, sine, cosine, ab, bc, ca, tx, ty, tz;
    float axis[3];
    float mag;

    axis[0] = ax;
    axis[1] = ay;
    axis[2] = az;
    mag = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (mag) {
        axis[0] /= mag;
        axis[1] /= mag;
        axis[2] /= mag;
    }

    radians = angle * myPi / 180.0f;
    sine = sin(radians);
    cosine = cos(radians);
    ab = axis[0] * axis[1] * (1 - cosine);
    bc = axis[1] * axis[2] * (1 - cosine);
    ca = axis[2] * axis[0] * (1 - cosine);
    tx = axis[0] * axis[0];
    ty = axis[1] * axis[1];
    tz = axis[2] * axis[2];

    m[0]  = tx + cosine * (1 - tx);
    m[1]  = ab + axis[2] * sine;
    m[2]  = ca - axis[1] * sine;
    m[3]  = 0.0f;
    m[4]  = ab - axis[2] * sine;
    m[5]  = ty + cosine * (1 - ty);
    m[6]  = bc + axis[0] * sine;
    m[7]  = 0.0f;
    m[8]  = ca + axis[1] * sine;
    m[9]  = bc - axis[0] * sine;
    m[10] = tz + cosine * (1 - tz);
    m[11] = 0;
    m[12] = 0;
    m[13] = 0;
    m[14] = 0;
    m[15] = 1;
}

void buildTranslateMatrix(float x, float y, float z, float m[16])
{
    m[0]  = 1;  m[1]  = 0;  m[2]  = 0;  m[3]  = x;
    m[4]  = 0;  m[5]  = 1;  m[6]  = 0;  m[7]  = y;
    m[8]  = 0;  m[9]  = 0;  m[10] = 1;  m[11] = z;
    m[12] = 0;  m[13] = 0;  m[14] = 0;  m[15] = 1;
}

void multMatrix(float dst[16],
                const float src1[16], const float src2[16])
{
    float tmp[16];
    int i, j;

    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            tmp[i*4+j] = src1[i*4+0] * src2[0*4+j] +
                src1[i*4+1] * src2[1*4+j] +
                src1[i*4+2] * src2[2*4+j] +
                src1[i*4+3] * src2[3*4+j];
        }
    }
    for (i=0; i<16; i++)
        dst[i] = tmp[i];
}

void transposeMatrix(float m[16])
{
    float tmp[16];
    int i,j;
    for(i=0; i<16; i++)
        tmp[i] = m[i];

    for(i=0; i<4; i++)
        for(j=0; j<4; j++)
            m[i*4+j] = tmp[i+j*4];
}


void multMatirxVecotr4(float dst[4], const float srcm[16], const float srcv[4])
{
	float t[4];
	t[0] = srcv[0];
	t[1] = srcv[1];
	t[2] = srcv[2];
	t[3] = srcv[3];

	dst[0] = srcm[0]*t[0] + srcm[1]*t[1] + srcm[2]*t[2] + srcm[3]*t[3];
	dst[1] = srcm[4]*t[0] + srcm[5]*t[1] + srcm[6]*t[2] + srcm[7]*t[3];
	dst[2] = srcm[8]*t[0] + srcm[9]*t[1] + srcm[10]*t[2] + srcm[11]*t[3];
	dst[3] = srcm[12]*t[0] + srcm[13]*t[1] + srcm[14]*t[2] + srcm[15]*t[3];
}

void multMatirxVecotr3(float dst[3], const float srcm[16], const float srcv[3])
{
	float t[3];
	t[0] = srcv[0];
	t[1] = srcv[1];
	t[2] = srcv[2];

	dst[0] = srcm[0] *t[0] + srcm[1] *t[1] + srcm[2] *t[2] + srcm[3];
	dst[1] = srcm[4] *t[0] + srcm[5] *t[1] + srcm[6] *t[2] + srcm[7];
	dst[2] = srcm[8] *t[0] + srcm[9] *t[1] + srcm[10]*t[2] + srcm[11];
}
