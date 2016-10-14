#pragma once

#include <algorithm>
#include <vector>
using namespace std;

class Element
{
public :
	int i;
	int j;
	double value;

	Element() { i = -1; j = -1; value = -1;}

	Element(int row, int col, double val)
	{
		i = row;
		j = col;
		value = val;
	}
};


class SparseMatrix
{
public:
	SparseMatrix(void);
	~SparseMatrix(void);
	SparseMatrix(int nRow, int nCol);

	void AddRow();

	vector<vector<Element>>::iterator Columns();

	int ColumnSize();

	vector<vector<Element>>::iterator Rows();

	int RowSize();


	Element AddElement(int i, int j, double value);
	Element AddElement(Element e);


	void Multiply(double* xIn, double* xOut);


	Element FindElement(int i, int j);

	void SortElement();
	Element AddValueTo(int i, int j, double value);
private:
	vector<vector<Element>> columns;
	vector<vector<Element>> rows;
	int m;
	int n;
};
