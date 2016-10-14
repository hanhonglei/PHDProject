#include "Miniball_dynamic_d.h"
// Class Implementations
// =====================

// Miniball_C
// --------


void Miniball_C::check_in (const Point& p)
{
	assert (d == p.dim());
	L.push_back(p);
}   

void Miniball_C::build ()
{
	B.reset();
	support_end = L.begin();
	pivot_mb (L.end());
}



void Miniball_C::mtf_mb (It i)
{
	support_end = L.begin();
	if ((B.size())==d+1) return;
	for (It k=L.begin(); k!=i;) {
		It j=k++;
		if (B.excess(*j) > 0) {
			if (B.push(*j)) {
				mtf_mb (j);
				B.pop();
				move_to_front(j);
			}
		}
	}
}


void Miniball_C::move_to_front (It j)
{
	if (support_end == j)
		support_end++;
	L.splice (L.begin(), L, j);
}



void Miniball_C::pivot_mb (It i)
{
	It t = ++L.begin();
	mtf_mb (t);
	double max_e, old_sqr_r = -1;
	do {
		It pivot;
		max_e = max_excess (t, i, pivot);
		if (max_e > 0) {
			t = support_end;
			if (t==pivot) ++t;
			old_sqr_r = B.squared_radius();
			B.push (*pivot);
			mtf_mb (support_end);
			B.pop();
			move_to_front (pivot);
		}
	} while ((max_e > 0) && (B.squared_radius() > old_sqr_r));
}



double Miniball_C::max_excess (It t, It i, It& pivot) const
{
	const double *c = B.center(), sqr_r = B.squared_radius();
	double e, max_e = 0;
	for (It k=t; k!=i; ++k) {
		const double *p = (*k).begin();
		e = -sqr_r;
		for (int j=0; j<d; ++j)
			e += mb_sqr(p[j]-c[j]);
		if (e > max_e) {
			max_e = e;
			pivot = k;
		}
	}
	return max_e;
}




Point Miniball_C::center () const
{
	return Point(d, B.center());
}


double Miniball_C::squared_radius () const
{
	return B.squared_radius();
}



int Miniball_C::nr_points () const
{
	return L.size();
}


Miniball_C::Cit Miniball_C::points_begin () const
{
	return L.begin();
}


Miniball_C::Cit Miniball_C::points_end () const
{
	return L.end();
} 

int Miniball_C::nr_support_points () const
{
	return B.support_size();
}


Miniball_C::Cit Miniball_C::support_points_begin () const
{
	return L.begin();
}


Miniball_C::Cit Miniball_C::support_points_end () const
{
	return support_end;
}

double Miniball_C::accuracy (double& slack) const
{
	double e, max_e = 0;
	int n_supp=0;
	Cit i;
	for (i=L.begin(); i!=support_end; ++i,++n_supp)
		if ((e = std::abs (B.excess (*i))) > max_e)
			max_e = e;

	// you've found a non-numerical problem if the following ever fails
	assert (n_supp == nr_support_points());

	for (i=support_end; i!=L.end(); ++i)
		if ((e = B.excess (*i)) > max_e)
			max_e = e;

	slack = B.slack();
	return (max_e/squared_radius());
}


bool Miniball_C::is_valid (double tolerance) const
{
	double slack;
	return ( (accuracy (slack) < tolerance) && (slack == 0) );
}   

// Miniball_C_b
// ----------


const double* Miniball_C_b::center () const
{
	return current_c;
}


double Miniball_C_b::squared_radius() const
{
	return current_sqr_r;
}


int Miniball_C_b::size() const
{
	return m;
}


int Miniball_C_b::support_size() const
{
	return s;
}


double Miniball_C_b::excess (const Point& p) const
{
	double e = -current_sqr_r;
	for (int k=0; k<d; ++k)
		e += mb_sqr(p[k]-current_c[k]);
	return e;
}




void Miniball_C_b::reset ()
{
	m = s = 0;
	// we misuse c[0] for the center of the empty sphere
	for (int j=0; j<d; ++j)
		c[0][j]=0;
	current_c = c[0];
	current_sqr_r = -1;
}


void Miniball_C_b::pop ()
{
	--m;
}



bool Miniball_C_b::push (const Point& p)
{
	int i, j;
	double eps = 1e-32;
	if (m==0) {
		for (i=0; i<d; ++i)
			q0[i] = p[i];
		for (i=0; i<d; ++i)
			c[0][i] = q0[i];
		sqr_r[0] = 0;
	} else {
		// set v_m to Q_m
		for (i=0; i<d; ++i)
			v[m][i] = p[i]-q0[i];

		// compute the a_{m,i}, i< m
		for (i=1; i<m; ++i) {
			a[m][i] = 0;
			for (j=0; j<d; ++j)
				a[m][i] += v[i][j] * v[m][j];
			a[m][i]*=(2/z[i]);
		}

		// update v_m to Q_m-\bar{Q}_m
		for (i=1; i<m; ++i) {
			for (j=0; j<d; ++j)
				v[m][j] -= a[m][i]*v[i][j];
		}

		// compute z_m
		z[m]=0;
		for (j=0; j<d; ++j)
			z[m] += mb_sqr(v[m][j]);
		z[m]*=2;

		// reject push if z_m too small
		if (z[m]<eps*current_sqr_r) {
			return false;
		}

		// update c, sqr_r
		double e = -sqr_r[m-1];
		for (i=0; i<d; ++i)
			e += mb_sqr(p[i]-c[m-1][i]);
		f[m]=e/z[m];

		for (i=0; i<d; ++i)
			c[m][i] = c[m-1][i]+f[m]*v[m][i];
		sqr_r[m] = sqr_r[m-1] + e*f[m]/2;
	}
	current_c = c[m];
	current_sqr_r = sqr_r[m];
	s = ++m;
	return true;
}


double Miniball_C_b::slack () const
{
	double* l = new double[d+1];
	double min_l=0;
	l[0] = 1;
	for (int i=s-1; i>0; --i) {
		l[i] = f[i];
		for (int k=s-1; k>i; --k)
			l[i]-=a[k][i]*l[k];
		if (l[i] < min_l) min_l = l[i];
		l[0] -= l[i];
	}
	if (l[0] < min_l) min_l = l[0];
	delete[] l;
	return ( (min_l < 0) ? -min_l : 0);
}

// Point
// -----

// Output

std::ostream& operator << (std::ostream& os, const Point& p)
{
	os << "(";
	int d = p.dim();
	for (int i=0; i<d-1; ++i)
		os << p[i] << ", ";
	os << p[d-1] << ")";
	return os;
}
