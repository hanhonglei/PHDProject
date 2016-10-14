#pragma once
#include "SparseMatrix.h"

class CCSMatrix
{
public:
	CCSMatrix(void);
	~CCSMatrix(void);

	CCSMatrix(SparseMatrix matrix);

	CCSMatrix(int _m, int _n);
	int* ColIndex();

	int ColumnSize();

	int NumNonZero();

	int* RowIndex();

	int RowSize();
	double* Values();

	int* rowIndex;
	double* values;
	int* colIndex;

	int nValuesNum;
	int m;
	int n;
};
