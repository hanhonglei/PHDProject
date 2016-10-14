#include "CSMatrix.h"

CCSMatrix::CCSMatrix(void)
{
	colIndex = NULL;
	rowIndex = NULL;
	values = NULL;
}

CCSMatrix::~CCSMatrix(void)
{
	if (colIndex!=NULL)
	{
		delete[] colIndex;
		colIndex = NULL;
	}

	if (rowIndex!=NULL)
	{
		delete[] rowIndex;
		rowIndex = NULL;
	}

	if (values!=NULL)
	{
		delete[] values;
		values = NULL;
	}
}


CCSMatrix::CCSMatrix(SparseMatrix matrix)
{
	colIndex = NULL;
	rowIndex = NULL;
	values = NULL;

	m = matrix.RowSize();
	n = matrix.ColumnSize();
	int nnz = 0;
	vector<vector<Element>>::iterator it = matrix.Columns();
	for (int i = 0; i < matrix.ColumnSize(); i++)
	{
		nnz += it[i].size();
	}
	//foreach (List<SparseMatrix.Element> col in matrix.Columns)
	//{
	//	nnz += col.Count;
	//}
	rowIndex = new int[nnz];
	colIndex = new int[n + 1];
	nValuesNum = nnz;
	values = new double[nValuesNum];
	int index = 0;
	int index2 = 0;
	colIndex[0] = 0;

	for (int i = 0; i < matrix.ColumnSize(); i++)
	{
		for (vector<Element>::iterator it2 = it[i].begin(); it2 != it[i].end(); it2++)
		{
			rowIndex[index] = it2->i;
			values[index] = it2->value;
			index++;
		}
		colIndex[++index2] = index;
	}
	//foreach (List<SparseMatrix.Element> col in matrix.Columns)
	//{
	//	foreach (SparseMatrix.Element e in col)
	//	{
	//		rowIndex[index] = e.i;
	//		values[index] = e.value;
	//		index++;
	//	}
	//	colIndex[++index2] = index;
	//}
}


CCSMatrix::CCSMatrix(int _m, int _n)
{
	colIndex = NULL;
	rowIndex = NULL;
	values = NULL;
	m = _m;
	n = _n;
	colIndex = new int[n + 1];
	colIndex[0] = 0;
}

int* CCSMatrix::ColIndex()
{
	return colIndex;
}

int CCSMatrix::ColumnSize()
{
	return n;
}

int CCSMatrix::NumNonZero()
{
	return nValuesNum/*values.Length*/;
}

int* CCSMatrix::RowIndex()
{
	return rowIndex;
}

int CCSMatrix::RowSize()
{
	return m;
}

double* CCSMatrix::Values()
{
	return values;
}
//////////////////////////////////////////////////////////////////////////
// ÔÝÊ±É¾µôµÄº¯Êý
//public void CG(double[] x, double[] b, double eps, int maxIter)
//{
//	double[] inv = new double[m];
//	double[] r = new double[m];
//	double[] d = new double[m];
//	double[] q = new double[m];
//	double[] s = new double[m];
//	double errNew = 0.0;
//	double err = 0.0;
//	double errOld = 0.0;
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
//		{
//			int row = rowIndex[j];
//			if (i == row)
//			{
//				inv[i] = 1.0 / values[j];
//			}
//		}
//	}
//	int iter = 0;
//	PreMultiply(x, r);
//	for (int i = 0; i < m; i++)
//	{
//		r[i] = b[i] - r[i];
//	}
//	for (int i = 0; i < m; i++)
//	{
//		d[i] = inv[i] * r[i];
//	}
//	for (int i = 0; i < m; i++)
//	{
//		errNew += r[i] * d[i];
//	}
//	err = errNew;
//	while ((iter < maxIter) && (errNew > (eps * err)))
//	{
//		PreMultiply(d, q);
//		double alpha = 0.0;
//		for (int i = 0; i < m; i++)
//		{
//			alpha += d[i] * q[i];
//		}
//		alpha = errNew / alpha;
//		for (int i = 0; i < m; i++)
//		{
//			x[i] += alpha * d[i];
//		}
//		if ((iter % 500) == 0)
//		{
//			PreMultiply(x, r);
//			for (int i = 0; i < m; i++)
//			{
//				r[i] = b[i] - r[i];
//			}
//		}
//		else
//		{
//			for (int i = 0; i < m; i++)
//			{
//				r[i] -= alpha * q[i];
//			}
//		}
//		for (int i = 0; i < m; i++)
//		{
//			s[i] = inv[i] * r[i];
//		}
//		errOld = errNew;
//		errNew = 0.0;
//		for (int i = 0; i < m; i++)
//		{
//			errNew += r[i] * s[i];
//		}
//		double beta = errNew / errOld;
//		for (int i = 0; i < m; i++)
//		{
//			d[i] = s[i] + (beta * d[i]);
//		}
//		iter++;
//	}
//}

//public bool Check(CCSMatrix B)
//{
//	CCSMatrix A = this;
//	if (A.rowIndex.Length != B.rowIndex.Length)
//	{
//		throw new Exception();
//	}
//	if (A.colIndex.Length != B.colIndex.Length)
//	{
//		throw new Exception();
//	}
//	if (A.values.Length != B.values.Length)
//	{
//		throw new Exception();
//	}
//	for (int i = 0; i < rowIndex.Length; i++)
//	{
//		if (A.rowIndex[i] != B.rowIndex[i])
//		{
//			throw new Exception();
//		}
//	}
//	for (int i = 0; i < colIndex.Length; i++)
//	{
//		if (A.colIndex[i] != B.colIndex[i])
//		{
//			throw new Exception();
//		}
//	}
//	for (int i = 0; i < values.Length; i++)
//	{
//		if (A.values[i] != B.values[i])
//		{
//			throw new Exception();
//		}
//	}
//	return true;
//}

//public CCSMatrix FastMultiply(CCSMatrix B)
//{
//	CCSMatrix A = this;
//	CCSMatrix C = new CCSMatrix(A.m, B.n);
//	Set<int> tmpIndex = new Set<int>();
//	double[] tmp = new double[A.m];
//	List<int> C_RowIndex = new List<int>();
//	List<double> C_Value = new List<double>();
//	for (int i = 0; i < A.m; i++)
//	{
//		tmp[i] = 0.0;
//	}
//	for (int j = 0; j < B.n; j++)
//	{
//		for (int col = B.colIndex[j]; col < B.colIndex[j + 1]; col++)
//		{
//			int k = B.rowIndex[col];
//			double valB = B.values[col];
//			if (k < A.ColumnSize)
//			{
//				for (int col2 = A.colIndex[k]; col2 < A.colIndex[k + 1]; col2++)
//				{
//					int k2 = A.rowIndex[col2];
//					double valA = A.values[col2];
//					tmpIndex.Add(k2);
//					tmp[k2] += valA * valB;
//				}
//			}
//		}
//		int[] t = tmpIndex.ToArray();
//		int count = 0;
//		Array.Sort<int>(t);
//		foreach (int k in t)
//		{
//			if (tmp[k] != 0.0)
//			{
//				C_RowIndex.Add(k);
//				C_Value.Add(tmp[k]);
//				tmp[k] = 0.0;
//				count++;
//			}
//		}
//		C.colIndex[j + 1] = C.colIndex[j] + count;
//		tmpIndex.Clear();
//	}
//	C.rowIndex = C_RowIndex.ToArray();
//	C_RowIndex = null;
//	C.values = C_Value.ToArray();
//	C_Value = null;
//	return C;
//}

/*	public void Multiply(double[] xIn, double[] xOut)
{
if ((xIn.Length < n) || (xOut.Length < m))
{
throw new ArgumentException();
}
for (int i = 0; i < m; i++)
{
xOut[i] = 0.0;
}
for (int i = 0; i < n; i++)
{
for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
{
int r = rowIndex[j];
xOut[r] += values[j] * xIn[i];
}
}
}*/

/*	public void PCG(double[] x, double[] b, double[] inv, double eps, int maxIter)
{
double[] r = new double[m];
double[] d = new double[m];
double[] q = new double[m];
double[] s = new double[m];
for (int i = 0; i < m; i++)
{
r[i] = b[i];
}
for (int i = 0; i < n; i++)
{
for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
{
r[rowIndex[j]] -= values[j] * x[i];
}
}
double new_err = 0.0;
for (int i = 0; i < m; i++)
{
d[i] = inv[i] * r[i];
new_err += d[i] * r[i];
}
double err = new_err;
for (int iter = 0; (iter < maxIter) && (new_err > ((eps * eps) * err)); iter++)
{
Multiply(d, q);
double tmp = 0.0;
for (int i = 0; i < m; i++)
{
tmp += d[i] * q[i];
}
double alpha = new_err / tmp;
for (int i = 0; i < m; i++)
{
x[i] += alpha * d[i];
}
for (int i = 0; i < m; i++)
{
if (double.IsNaN(x[i]))
{
throw new Exception();
}
}
if ((iter % 50) == 0)
{
for (int i = 0; i < m; i++)
{
r[i] = b[i];
}
for (int i = 0; i < n; i++)
{
for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
{
r[rowIndex[j]] -= values[j] * x[i];
}
}
}
else
{
for (int i = 0; i < m; i++)
{
r[i] -= alpha * q[i];
}
}
for (int i = 0; i < m; i++)
{
s[i] = inv[i] * r[i];
}
double old_err = new_err;
new_err = 0.0;
for (int i = 0; i < m; i++)
{
new_err += r[i] * s[i];
}
double beta = new_err / old_err;
for (int i = 0; i < m; i++)
{
d[i] = s[i] + (beta * d[i]);
}
}
}*/

//public void PreMultiply(double[] xIn, double[] xOut)
//{
//	if ((xIn.Length < m) || (xOut.Length < n))
//	{
//		throw new ArgumentException();
//	}
//	for (int i = 0; i < n; i++)
//	{
//		xOut[i] = 0.0;
//	}
//	for (int i = 0; i < n; i++)
//	{
//		double sum = 0.0;
//		for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
//		{
//			int r = rowIndex[j];
//			sum += values[j] * xIn[r];
//		}
//		xOut[i] = sum;
//	}
//}

//public void PreMultiply(double[] xIn, double[] xOut, int[] index)
//{
//	if ((xIn.Length < m) || (xOut.Length < n))
//	{
//		throw new ArgumentException();
//	}
//	foreach (int i in index)
//	{
//		xOut[i] = 0.0;
//	}
//	foreach (int i in index)
//	{
//		double sum = 0.0;
//		for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
//		{
//			int r = rowIndex[j];
//			sum += values[j] * xIn[r];
//		}
//		xOut[i] = sum;
//	}
//}

//public void PreMultiply(double[] xIn, double[] xOut, int startIndex, bool init)
//{
//	if ((xIn.Length < m) || (xOut.Length < n))
//	{
//		throw new ArgumentException();
//	}
//	if (init)
//	{
//		for (int i = 0; i < n; i++)
//		{
//			xOut[i] = 0.0;
//		}
//	}
//	for (int i = 0; i < n; i++)
//	{
//		double sum = 0.0;
//		for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
//		{
//			int r = rowIndex[j];
//			sum += values[j] * xIn[r + startIndex];
//		}
//		xOut[i] += sum;
//	}
//}

//public void PreMultiply(double[] xIn, double[] xOut, int inStart, int outStart, bool init)
//{
//	if ((xIn.Length < m) || (xOut.Length < n))
//	{
//		throw new ArgumentException();
//	}
//	if (init)
//	{
//		for (int i = 0; i < n; i++)
//		{
//			xOut[i + outStart] = 0.0;
//		}
//	}
//	for (int i = 0; i < n; i++)
//	{
//		double sum = 0.0;
//		for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
//		{
//			int r = rowIndex[j];
//			sum += values[j] * xIn[r + inStart];
//		}
//		xOut[i + outStart] += sum;
//	}
//}

//public void PreMultiplyOffset(double[] xIn, double[] xOut, int startOut, int offsetOut)
//{
//	int i = 0;
//	for (i = startOut; i < (n + offsetOut); i += offsetOut)
//	{
//		xOut[i] = 0.0;
//	}
//	i = 0;
//	for (int k = startOut; i < n; k += offsetOut)
//	{
//		double sum = 0.0;
//		for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
//		{
//			int r = rowIndex[j];
//			sum += values[j] * xIn[r];
//		}
//		xOut[k] = sum;
//		i++;
//	}
//}

/*	public CCSMatrix Transpose()
{
CCSMatrix C = new CCSMatrix(n, m);
int[] rowCount = new int[m];
for (int i = 0; i < rowCount.Length; i++)
{
rowCount[i] = 0;
}
for (int i = 0; i < rowIndex.Length; i++)
{
rowCount[rowIndex[i]]++;
}
C.ColIndex[0] = 0;
for (int i = 0; i < m; i++)
{
C.ColIndex[i + 1] = C.ColIndex[i] + rowCount[i];
rowCount[i] = C.ColIndex[i];
}
C.values = new double[NumNonZero];
C.rowIndex = new int[NumNonZero];
for (int i = 0; i < n; i++)
{
for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
{
int k = rowIndex[j];
C.values[rowCount[k]] = values[j];
C.rowIndex[rowCount[k]] = i;
rowCount[k]++;
}
}
return C;
}*/

/*	public void Write(StreamWriter sw)
{
for (int i = 0; i < n; i++)
{
for (int j = colIndex[i]; j < colIndex[i + 1]; j++)
{
sw.WriteLine(string.Concat(new object[] { rowIndex[j] + 1, " ", i + 1, " ", values[j].ToString() }));
if (i != rowIndex[j])
{
sw.WriteLine(string.Concat(new object[] { i + 1, " ", rowIndex[j] + 1, " ", values[j].ToString() }));
}
}
}
}*/
//
//CCSMatrix::CCSMatrix(double[,] matrix)
//{
//	m = matrix.GetLength(0);
//	n = matrix.GetLength(1);
//	int nnz = 0;
//	for (int i = 0; i < m; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			if (matrix[i, j] != 0.0)
//			{
//				nnz++;
//			}
//		}
//	}
//	rowIndex = new int[nnz];
//	colIndex = new int[n + 1];
//	nValuesNum = nnz;
//	values = new double[nValuesNum];
//	int index = 0;
//	int index2 = 0;
//	colIndex[0] = 0;
//	for (int j = 0; j < n; j++)
//	{
//		for (int i = 0; i < n; i++)
//		{
//			if (matrix[i, j] != 0.0)
//			{
//				rowIndex[index] = i;
//				values[index] = matrix[i, j];
//				index++;
//			}
//		}
//		colIndex[++index2] = index;
//	}
//}

//
//public CCSMatrix(SparseMatrix matrix, bool transponse)
//{
//	m = matrix.ColumnSize;
//	n = matrix.RowSize;
//	int nnz = 0;
//	foreach (List<SparseMatrix.Element> col in matrix.Columns)
//	{
//		nnz += col.Count;
//	}
//	rowIndex = new int[nnz];
//	colIndex = new int[n + 1];
//	nValuesNum = nnz;
//	values = new double[nValuesNum];
//	int index = 0;
//	int index2 = 0;
//	colIndex[0] = 0;
//	foreach (List<SparseMatrix.Element> row in matrix.Rows)
//	{
//		foreach (SparseMatrix.Element e in row)
//		{
//			rowIndex[index] = e.j;
//			values[index] = e.value;
//			index++;
//		}
//		colIndex[++index2] = index;
//	}
//}
