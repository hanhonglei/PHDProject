//    Copright (C) 1999-2006, Bernd Gaertner
//    $Revision: 1.4 $
//    $Date: 2008/02/11 11:50:05 $
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA,
//    or download the License terms from prep.ai.mit.edu/pub/gnu/COPYING-2.0.
//
//    Contact:
//    --------
//    Bernd Gaertner
//    Institute of Theoretical Computer Science 
//    ETH Zuerich
//    CAB G32.2
//    CH-8092 Zuerich, Switzerland
//    http://www.inf.ethz.ch/personal/gaertner

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <iostream>
#include <list>

// Functions
// =========
inline double mb_sqr (double r) {return r*r;}
 
// Class Declarations
// ==================

// smallest enclosing ball of a set of n points in dimension d
class Miniball_C;

// smallest ball with a set of n <= d+1 points *on the boundary*
class Miniball_C_b; 

// point in dimension d
class Point;


// Class Definitions
// =================

// Miniball_C_b
// ----------
class Miniball_C_b {
private:
    
  // data members
  int                 d;      // dimension
  int                 m, s;   // size and number of support points
  double*             q0;
    
  double*             z;
  double*             f;
  double**            v;
  double**            a;
    
  double**            c;
  double*             sqr_r;
    
  double*             current_c;      // refers to some c[j]
  double              current_sqr_r;  

  // if you want the copy constructor and assignment operator, please
  // properly implement them yourself. The default copy/assignment
  // semantics fails since a Miniball_C_b object stores pointers to
  // dynamically allocated memory
  Miniball_C_b (const Miniball_C_b& mb);
  Miniball_C_b& operator=(const Miniball_C_b& mb);  

public:
  Miniball_C_b(int dim) 
    : d(dim) 
  {
    q0 = new double[d];
    z = new double[d+1];
    f = new double[d+1];
    v = new double*[d+1];
    for (int i=0; i<d+1; ++i) v[i] =  new double[d];
    a = new double*[d+1];
    for (int i=0; i<d+1; ++i) a[i] =  new double[d];   
    c = new double*[d+1];
    for (int i=0; i<d+1; ++i) c[i] =  new double[d];
    sqr_r = new double[d+1];
    reset();
  }

  ~Miniball_C_b() 
  {
    delete[] sqr_r;
    for (int i=0; i<d+1; ++i) delete[] c[i];
    delete[] c;
    for (int i=0; i<d+1; ++i) delete[] a[i];
    delete[] a;
    for (int i=0; i<d+1; ++i) delete[] v[i];
    delete[] v;
    delete[] f;
    delete[] z;
    delete[] q0;
  }
    
  // access
  const double*       center() const;
  double              squared_radius() const;
  int                 size() const;
  int                 support_size() const;
  double              excess (const Point& p) const;
    
  // modification
  void                reset(); // generates empty sphere with m=s=0
  bool                push (const Point& p);
  void                pop ();
    
  // checking
  double              slack() const;
};

// Miniball_C
// --------
class Miniball_C {
public:
  // types
  typedef std::list<Point>::iterator         It;
  typedef std::list<Point>::const_iterator   Cit;
    
private:
  // data members
  int              d;            // dimension
  std::list<Point> L;            // internal point set
  Miniball_C_b       B;            // the current ball
  It               support_end;  // past-the-end iterator of support set
    
  // private methods
  void        mtf_mb (It k);
  void        pivot_mb (It k);
  void        move_to_front (It j);
  double      max_excess (It t, It i, It& pivot) const;  

public:
  // creates an empty ball
  Miniball_C(int dim) : d(dim), B(dim) {}

  // copies p to the internal point set
  void        check_in (const Point& p);

  // builds the smallest enclosing ball of the internal point set
  void        build ();
    
  // returns center of the ball (undefined if ball is empty)
  Point       center() const;

  // returns squared_radius of the ball (-1 if ball is empty)
  double      squared_radius () const;

  // returns size of internal point set
  int         nr_points () const;

  // returns begin- and past-the-end iterators for internal point set
  Cit         points_begin () const;
  Cit         points_end () const;

  // returns size of support point set; this set has the following properties:
  // - there are at most d+1 support points, 
  // - all support points are on the boundary of the computed ball, and
  // - the smallest enclosing ball of the support point set equals the 
  //   smallest enclosing ball of the internal point set
  int         nr_support_points () const;

  // returns begin- and past-the-end iterators for internal point set
  Cit         support_points_begin () const;
  Cit         support_points_end () const;
    
  // assesses the quality of the computed ball. The return value is the
  // maximum squared distance of any support point or point outside the 
  // ball to the boundary of the ball, divided by the squared radius of
  // the ball. If everything went fine, this will be less than e-15 and
  // says that the computed ball approximately contains all the internal
  // points and has all the support points on the boundary.
  // 
  // The slack parameter that is set by the method says something about
  // whether the computed ball is really the *smallest* enclosing ball 
  // of the support points; if everything went fine, this value will be 0; 
  // a positive value may indicate that the ball is not smallest possible,
  // with the deviation from optimality growing with the slack
  double      accuracy (double& slack) const;

  // returns true if the accuracy is below the given tolerance and the
  // slack is 0
  bool        is_valid (double tolerance = 1e-15) const;
};
    
// Point (inline)
// --------------

class Point {
private:
  int d; 
  double* coord;
   
public:
  // default
  Point(int dim)
    : d (dim), coord(new double[dim])
  {}

  ~Point ()
  {
    delete[] coord;
  }
   
  // copy from Point
  Point (const Point& p)
    : d (p.dim()), coord(new double[p.dim()])
  {
    for (int i=0; i<d; ++i)
      coord[i] = p.coord[i];
  }
   
  // copy from double*
  Point (int dim, const double* p)
    : d (dim), coord(new double[dim])
  {
    for (int i=0; i<d; ++i)
      coord[i] = p[i];
  }
   
  // assignment
  Point& operator = (const Point& p)
  {
    assert (d == p.dim());
    if (this != &p)
      for (int i=0; i<d; ++i)
	coord[i] = p.coord[i];
    return *this;
  }
   
  // coordinate access
  double& operator [] (int i)
  {
    return coord[i];
  }
  const double& operator [] (int i) const
  {
    return coord[i];
  }
  const double* begin() const
  {
    return coord;
  }
  const double* end() const
  {
    return coord+d;
  }

  // dimension access
  int dim() const
  {
    return d;
  }
};
   
