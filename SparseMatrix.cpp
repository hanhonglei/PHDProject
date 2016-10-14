#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(void)
{
}

SparseMatrix::~SparseMatrix(void)
{
}


bool ColumnComparer(Element e1, Element e2)
{
	return e1.i < e2.i;
}


bool RowComparer(Element e1, Element e2)
{
	return e1.j < e2.j;
}


void SparseMatrix::AddRow()
{
	//rows.Add(new List<Element>());
	rows.resize(rows.size()+1);
	m++;
}

Element SparseMatrix::AddValueTo(int i, int j, double value)
{
	Element e = FindElement(i, j);
	if (e.i == e.j && e.i == -1)
	{
		e = Element(i, j, value);
		AddElement(e);
	}
	else
	{
		for (vector<Element>::iterator it=rows[i].begin(); it != rows[i].end();++it)
		{
			if (it->j == j)
			{
				it->value += value;
				break;
			}
		}
		e.value += value;
	}
	return e;
}


vector<vector<Element>>::iterator SparseMatrix::Columns()
{

	return columns.begin();
}

int SparseMatrix::ColumnSize()
{
	return n;
}

vector<vector<Element>>::iterator SparseMatrix::Rows()
{
	return rows.begin();
}

int SparseMatrix::RowSize()
{
	return m;
}


SparseMatrix::SparseMatrix(int nRow, int nCol)
{
	m = nRow;
	n = nCol;
	//rows = new List<List<Element>>(m);
	//columns = new List<List<Element>>(n);
	rows.resize(nRow);
	columns.resize(nCol);
	//for (int i = 0; i < m; i++)
	//{
		rows.resize(m);
		//rows.Add(new List<Element>());
	//}
	//for (int i = 0; i < n; i++)
	//{
		columns.resize(n);
		//columns.Add(new List<Element>());
	//}
}


Element SparseMatrix::AddElement(int i, int j, double value)
{

	//List<Element> r = rows[i];
	//List<Element> c = columns[j];
	Element e = Element(i, j, value);
	rows[i].push_back(e);
	columns[j].push_back(e);
	//r.Add(e);
	//c.Add(e);
	return e;
}


void SparseMatrix::Multiply(double* xIn, double* xOut)
{
	//if ((xIn.Length < n) || (xOut.Length < m))
	//{
	//	throw new ArgumentException();
	//}
	for (int i = 0; i < m; i++)
	{
		//List<Element> r = rows[i];
		double sum = 0.0;

		//foreach (Element e in r)
		//{
		//	sum += e.value * xIn[e.j];
		//}
		for (vector<Element>::iterator it =rows[i].begin() ; it < rows[i].end(); it++)
		{
			sum += it->value * xIn[it->j];

		}
		xOut[i] = sum;
	}
}


Element SparseMatrix::FindElement(int i, int j)
{
	//List<Element> rr = rows[i];
	//foreach (Element e in rr)
	//{
	//	if (e.j == j)
	//	{
	//		return e;
	//	}
	//}
	Element e;
	for (vector<Element>::iterator it=rows[i].begin(); it != rows[i].end();++it)
	{
		if (it->j == j)
		{
			return *it;
		}
	}
	return e;
}

void SparseMatrix::SortElement()
{
	//RowComparer rComparer = new RowComparer();
	//ColumnComparer cComparer = new ColumnComparer();
	//foreach (List<Element> r in rows)
	//{
	//	r.Sort(rComparer);
	//}
	//foreach (List<Element> c in columns)
	//{
	//	c.Sort(cComparer);
	//}
	for (vector<vector<Element>>::iterator it=rows.begin(); it != rows.end();++it)
	{
		sort(it->begin(), it->end(), RowComparer);
	}
	for (vector<vector<Element>>::iterator it=columns.begin(); it != columns.end();++it)
	{
		sort(it->begin(), it->end(), ColumnComparer);
	}
}

Element SparseMatrix::AddElement(Element e)
{
	rows[e.i].push_back(e);
	columns[e.j].push_back(e);
	//List<Element> r = rows[e.i];
	//List<Element> c = columns[e.j];
	//r.Add(e);
	//c.Add(e);
	return e;
}

/*
public SparseMatrix(SparseMatrix right)
{
	m = right.m;
	n = right.n;
	rows = new List<List<Element>>(m);
	columns = new List<List<Element>>(n);
	for (int i = 0; i < m; i++)
	{
		rows.Add(new List<Element>());
	}
	for (int i = 0; i < n; i++)
	{
		columns.Add(new List<Element>());
	}
	foreach (List<Element> vector in right.Rows)
	{
		foreach (Element e in vector)
		{
			AddElement(e.i, e.j, e.value);
		}
	}
}

public SparseMatrix(int m, int n, int nElements)
{
	m = m;
	n = n;
	rows = new List<List<Element>>(m);
	columns = new List<List<Element>>(n);
	for (int i = 0; i < m; i++)
	{
		rows.Add(new List<Element>(nElements));
	}
	for (int i = 0; i < n; i++)
	{
		columns.Add(new List<Element>(nElements));
	}
}

public SparseMatrix Add(SparseMatrix right)
{
	if ((m != right.m) || (n != right.n))
	{
		throw new ArgumentException();
	}
	SparseMatrix ret = new SparseMatrix(m, m);
	for (int i = 0; i < m; i++)
	{
		List<Element> r1 = Rows[i];
		List<Element> r2 = right.Rows[i];
		int c1 = 0;
		int c2 = 0;
		while ((c1 < r1.Count) && (c2 < r2.Count))
		{
			Element e1 = r1[c1];
			Element e2 = r2[c2];
			if (e1.j < e2.j)
			{
				c1++;
				ret.AddElement(i, e1.j, e1.value);
			}
			else
			{
				if (e1.j > e2.j)
				{
					c2++;
					ret.AddElement(i, e2.j, e2.value);
					continue;
				}
				ret.AddElement(i, e1.j, e1.value + e2.value);
				c1++;
				c2++;
			}
		}
		while (c1 < r1.Count)
		{
			Element e = r1[c1];
			ret.AddElement(e.i, e.j, e.value);
			c1++;
		}
		while (c2 < r2.Count)
		{
			Element e = r2[c2];
			ret.AddElement(e.i, e.j, e.value);
			c2++;
		}
	}
	return ret;
}

private static void Add(double[] w, double[] u, double[] v)
{
	if ((u.Length != v.Length) || (v.Length != w.Length))
	{
		throw new ArgumentException();
	}
	for (int i = 0; i < u.Length; i++)
	{
		w[i] = u[i] + v[i];
	}
}

public void AddColumn()
{
	columns.Add(new List<Element>());
	n++;
}

public Element AddElement(Element e)
{
	List<Element> r = rows[e.i];
	List<Element> c = columns[e.j];
	r.Add(e);
	c.Add(e);
	return e;
}

public Element AddElementIfNotExist(int i, int j, double value)
{
	if (FindElement(i, j) == null)
	{
		return AddElement(i, j, value);
	}
	return null;
}

private static void Assign(double[] u, double[] v)
{
	if (u.Length != v.Length)
	{
		throw new ArgumentException();
	}
	for (int i = 0; i < u.Length; i++)
	{
		u[i] = v[i];
	}
}

public bool CheckElements()
{
	foreach (List<Element> r in rows)
	{
		foreach (Element e in r)
		{
			if ((double.IsInfinity(e.value) || double.IsNaN(e.value)) || (e.value == 0.0))
			{
				return false;
			}
		}
	}
	return true;
}

public SparseMatrix ConcatRows(SparseMatrix right)
{
	if (ColumnSize != right.ColumnSize)
	{
		throw new ArgumentException();
	}
	SparseMatrix m = new SparseMatrix(RowSize + right.RowSize, ColumnSize);
	foreach (List<Element> r in rows)
	{
		foreach (Element e in r)
		{
			m.AddElement(e.i, e.j, e.value);
		}
	}
	int r_base = RowSize;
	foreach (List<Element> r in right.rows)
	{
		foreach (Element e in r)
		{
			m.AddElement(r_base + e.i, e.j, e.value);
		}
	}
	return m;
}

public static void ConjugateGradientsMethod(SparseMatrix A, double[] b, double[] x, int iter, double tolerance)
{
	int n = A.ColumnSize;
	int rn = (int) Math.Sqrt((double) n);
	if (A.RowSize != A.ColumnSize)
	{
		throw new ArgumentException();
	}
	if ((b.Length != n) || (x.Length != n))
	{
		throw new ArgumentException();
	}
	double[] r = new double[n];
	double[] d = new double[n];
	double[] q = new double[n];
	double[] t = new double[n];
	int i = 0;
	A.Multiply(x, r);
	Subtract(r, b, r);
	Assign(d, r);
	double newError = Dot(r, r);
	while ((i < iter) && (newError > tolerance))
	{
		A.Multiply(d, q);
		double alpha = newError / Dot(d, q);
		Scale(t, d, alpha);
		Add(x, x, t);
		if ((i % rn) == 0)
		{
			A.Multiply(x, r);
			Subtract(r, b, r);
		}
		else
		{
			Scale(t, q, alpha);
			Subtract(r, r, t);
		}
		double oldError = newError;
		newError = Dot(r, r);
		if (newError < tolerance)
		{
			A.Multiply(x, r);
			Subtract(r, b, r);
			newError = Dot(r, r);
		}
		double beta = newError / oldError;
		Scale(d, d, beta);
		Add(d, d, r);
		i++;
	}
}

public static void ConjugateGradientsMethod2(SparseMatrix A, double[] b, double[] x, int iter, double tolerance)
{
	int n = A.ColumnSize;
	int rn = (int) Math.Sqrt((double) n);
	if (A.RowSize != A.ColumnSize)
	{
		throw new ArgumentException();
	}
	if ((b.Length != n) || (x.Length != n))
	{
		throw new ArgumentException();
	}
	double[] r1 = new double[n];
	double[] r2 = new double[n];
	double[] d1 = new double[n];
	double[] d2 = new double[n];
	double[] q1 = new double[n];
	double[] q2 = new double[n];
	double[] t = new double[n];
	int i = 0;
	A.Multiply(x, r1);
	Subtract(r1, b, r1);
	Assign(r2, r1);
	Assign(d1, r1);
	Assign(d2, r2);
	double newError = Dot(r2, r1);
	while ((i < iter) && (newError > tolerance))
	{
		A.Multiply(d1, q1);
		A.PreMultiply(d2, q2);
		double alpha = newError / Dot(d2, q1);
		Scale(t, d1, alpha);
		Add(x, x, t);
		if ((i % rn) == 0)
		{
			A.Multiply(x, r1);
			Subtract(r1, b, r1);
			Assign(r2, r1);
		}
		else
		{
			Scale(t, q1, alpha);
			Subtract(r1, r1, t);
			Scale(t, q2, alpha);
			Subtract(r2, r2, t);
		}
		double oldError = newError;
		newError = Dot(r2, r1);
		if (newError < tolerance)
		{
			A.Multiply(x, r1);
			Subtract(r1, b, r1);
			Assign(r2, r1);
			newError = Dot(r2, r1);
		}
		double beta = newError / oldError;
		Scale(d1, d1, beta);
		Add(d1, d1, r1);
		Scale(d2, d2, beta);
		Add(d2, d2, r2);
		i++;
	}
}

public static void ConjugateGradientsMethod3(SparseMatrix A, double[] b, double[] x, bool[] boundary, int iter, double tolerance)
{
	int n = A.ColumnSize;
	int rn = (int) Math.Sqrt((double) n);
	if (A.RowSize != A.ColumnSize)
	{
		throw new ArgumentException();
	}
	if ((b.Length != n) || (x.Length != n))
	{
		throw new ArgumentException();
	}
	double[] r = new double[n];
	double[] d = new double[n];
	double[] q = new double[n];
	double[] t = new double[n];
	int i = 0;
	A.Multiply(x, r);
	Subtract(r, b, r);
	Assign(d, r);
	double newError = Dot(r, r);
	while ((i < iter) && (newError > tolerance))
	{
		A.Multiply(d, q);
		double alpha = newError / Dot(d, q);
		Scale(t, d, alpha);
		Add(x, x, t);
		if ((i % rn) == 0)
		{
			A.Multiply(x, t);
			Subtract(r, b, t);
		}
		else
		{
			Scale(t, q, alpha);
			Subtract(r, r, t);
		}
		double oldError = newError;
		newError = Dot(r, r);
		if (newError < tolerance)
		{
			A.Multiply(x, r);
			Subtract(r, b, r);
			newError = Dot(r, r);
		}
		double beta = newError / oldError;
		Scale(d, d, beta);
		Add(d, d, r);
		i++;
	}
}

private static double Dot(double[] u, double[] v)
{
	if (u.Length != v.Length)
	{
		throw new ArgumentException();
	}
	double sum = 0.0;
	for (int i = 0; i < u.Length; i++)
	{
		sum += u[i] * v[i];
	}
	return sum;
}

public override bool Equals(object obj)
{
	SparseMatrix right = obj as SparseMatrix;
	if (obj == null)
	{
		return false;
	}
	if (right.m != m)
	{
		return false;
	}
	if (right.n != n)
	{
		return false;
	}
	for (int i = 0; i < n; i++)
	{
		List<Element> c1 = columns[i];
		List<Element> c2 = right.columns[i];
		if (c1.Count != c2.Count)
		{
			return false;
		}
		for (int j = 0; j < c1.Count; j++)
		{
			Element e1 = c1[j];
			Element e2 = c2[j];
			if (e1.j != e2.j)
			{
				return false;
			}
			if (e1.value != e2.value)
			{
				return false;
			}
		}
	}
	return true;
}


public List<Element> GetColumn(int index)
{
	return columns[index];
}

public int[][] GetColumnIndex()
{
	int[][] arr = new int[n][];
	for (int i = 0; i < n; i++)
	{
		arr[i] = new int[columns[i].Count];
		int j = 0;
		foreach (Element e in columns[i])
		{
			arr[i][j++] = e.i;
		}
	}
	return arr;
}

public double[] GetDiagonalPreconditionor()
{
	if (m != n)
	{
		return null;
	}
	double[] ret = new double[n];
	for (int i = 0; i < n; i++)
	{
		Element d = FindElement(i, i);
		if (d == null)
		{
			ret[i] = 1.0;
		}
		else if (d.value == 0.0)
		{
			ret[i] = 1.0;
		}
		else
		{
			ret[i] = 1.0 / d.value;
		}
	}
	return ret;
}

public override int GetHashCode()
{
	return (base.GetHashCode() + NumOfElements());
}

public List<Element> GetRow(int index)
{
	return rows[index];
}

public int[][] GetRowIndex()
{
	int[][] arr = new int[m][];
	for (int i = 0; i < m; i++)
	{
		arr[i] = new int[rows[i].Count];
		int j = 0;
		foreach (Element e in rows[i])
		{
			arr[i][j++] = e.j;
		}
	}
	return arr;
}

public bool IsSymmetric()
{
	if (m != n)
	{
		return false;
	}
	for (int i = 0; i < m; i++)
	{
		List<Element> row = GetRow(i);
		List<Element> col = GetColumn(i);
		if (row.Count != col.Count)
		{
			return false;
		}
		for (int j = 0; j < row.Count; j++)
		{
			Element e1 = row[j];
			Element e2 = col[j];
			if (e1.i != e2.j)
			{
				return false;
			}
			if (e1.j != e2.i)
			{
				return false;
			}
		}
	}
	return true;
}

public static void JacobiMethod(SparseMatrix A, double[] b, double[] x, int iter, double tolerance)
{
	int n = A.ColumnSize;
	if (A.ColumnSize != A.RowSize)
	{
		throw new ArgumentException();
	}
	if ((b.Length != n) || (x.Length != n))
	{
		throw new ArgumentException();
	}
	double[] r = new double[n];
	double[] d = new double[n];
	double[] t = new double[n];
	for (int i = 0; i < n; i++)
	{
		Element e = A.FindElement(i, i);
		if (e != null)
		{
			d[i] = 0.75 / e.value;
		}
		else
		{
			d[i] = 0.0;
		}
	}
	A.Multiply(x, t);
	Subtract(r, b, t);
	double error = Dot(r, r);
	for (int count = 0; (count < iter) && (error > tolerance); count++)
	{
		for (int i = 0; i < n; i++)
		{
			t[i] = r[i] * d[i];
		}
		Add(x, x, t);
		A.Multiply(x, t);
		Subtract(r, b, t);
		error = Dot(r, r);
	}
}

public SparseMatrix Multiply(SparseMatrix right)
{
	if (n != right.m)
	{
		throw new ArgumentException();
	}
	SparseMatrix ret = new SparseMatrix(m, right.n);
	for (int i = 0; i < Rows.Count; i++)
	{
		List<Element> rr = Rows[i];
		for (int j = 0; j < right.Columns.Count; j++)
		{
			List<Element> cc = right.Columns[j];
			int c1 = 0;
			int c2 = 0;
			double sum = 0.0;
			bool used = false;
			while ((c1 < rr.Count) && (c2 < cc.Count))
			{
				Element e1 = rr[c1];
				Element e2 = cc[c2];
				if (e1.j < e2.i)
				{
					c1++;
				}
				else
				{
					if (e1.j > e2.i)
					{
						c2++;
						continue;
					}
					sum += e1.value * e2.value;
					c1++;
					c2++;
					used = true;
				}
			}
			if (used)
			{
				ret.AddElement(i, j, sum);
			}
		}
	}
	return ret;
}

public void Multiply(double[] xIn, int indexIn, double[] xOut, int indexOut)
{
	if (((xIn.Length - indexIn) < n) || ((xOut.Length - indexOut) < m))
	{
		throw new ArgumentException();
	}
	for (int i = 0; i < m; i++)
	{
		List<Element> r = rows[i];
		double sum = 0.0;
		foreach (Element e in r)
		{
			sum += e.value * xIn[e.j + indexIn];
		}
		xOut[i + indexOut] = sum;
	}
}

public int NumOfElements()
{
	int count = 0;
	if (m < n)
	{
		foreach (List<Element> r in rows)
		{
			count += r.Count;
		}
		return count;
	}
	foreach (List<Element> c in columns)
	{
		count += c.Count;
	}
	return count;
}

public void PreMultiply(double[] xIn, double[] xOut)
{
	if ((xIn.Length < m) || (xOut.Length < n))
	{
		throw new ArgumentException();
	}
	for (int j = 0; j < n; j++)
	{
		List<Element> c = columns[j];
		double sum = 0.0;
		foreach (Element e in c)
		{
			sum += e.value * xIn[e.i];
		}
		xOut[j] = sum;
	}
}

public void Scale(double s)
{
	foreach (List<Element> vector in rows)
	{
		foreach (Element e in vector)
		{
			e.value *= s;
		}
	}
}

private static void Scale(double[] w, double[] u, double s)
{
	if (u.Length != w.Length)
	{
		throw new ArgumentException();
	}
	for (int i = 0; i < u.Length; i++)
	{
		w[i] = u[i] * s;
	}
}


private static void Subtract(double[] w, double[] u, double[] v)
{
	if ((u.Length != v.Length) || (v.Length != w.Length))
	{
		throw new ArgumentException();
	}
	for (int i = 0; i < u.Length; i++)
	{
		w[i] = u[i] - v[i];
	}
}

public SparseMatrix Transpose()
{
	SparseMatrix ret = new SparseMatrix(this);
	int t = ret.m;
	ret.m = ret.n;
	ret.n = t;
	List<List<Element>> tmp = ret.rows;
	ret.rows = ret.columns;
	ret.columns = tmp;
	foreach (List<Element> r in ret.rows)
	{
		foreach (Element e in r)
		{
			t = e.i;
			e.i = e.j;
			e.j = t;
		}
	}
	return ret;
}

*/

