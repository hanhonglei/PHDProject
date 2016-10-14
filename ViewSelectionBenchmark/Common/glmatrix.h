/************************************************************************/
/*          These functions are from Nv's Cg Examples,                  */
/*          For Building and mul-matrix for Cg Usage.                   */
/************************************************************************/

//All this matrix are in row-major mode, for Cg usage

void buildIdentityMatrix(float m[16]);

void buildOrthoMatrix(double left, double right,
                      double bottom, double top,
                      double near, double far,
                      float m[16]);

void buildPerspectiveMatrix(double fieldOfView,
                            double aspectRatio,
                            double zNear, double zFar,
                            float m[16]);

void buildLookAtMatrix(double eyex, double eyey, double eyez,
                       double centerx, double centery, double centerz,
                       double upx, double upy, double upz,
                       float m[16]);

void buildRotateMatrix(float angle,
                      float ax, float ay, float az,
                      float m[16]);

void buildTranslateMatrix(float x, float y, float z, float m[16]);

void multMatrix(float dst[16],
                const float src1[16], const float src2[16]);

void transposeMatrix(float m[16]);
void multMatirxVecotr4(float dst[4], const float srcm[16], const float srcv[4]);
void multMatirxVecotr3(float dst[3], const float srcm[16], const float srcv[3]);