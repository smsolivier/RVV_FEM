#pragma once 

#include "General.hpp"

namespace fem 
{

/// represent a location in DIM dimensions 
class Point {
public:
	/// initialize to zero 
	Point() {for (int d=0; d<DIM; d++) _x[d] = 0.; }
	/// initializer list constructor 
	Point(std::initializer_list<double> list); 
	/// index into dth dimension 
	double& operator[](int d) {
		CHECK(d < DIM); 
		return _x[d]; 
	}
	/// const index into dth dimension
	double operator[](int d) const {
		CHECK(d < DIM); 
		return _x[d]; 
	}
	/// add two points together 
	void operator+=(const Point& x) {
		for (int i=0; i<DIM; i++) {
			(*this)[i] += x[i]; 
		}
	}

	/// divide all components of this by a double 
	void operator/=(double x) {
		for (int i=0; i<DIM; i++) {
			(*this)[i] /= x; 
		}
	}

	/// subtract two points 
	Point operator-(const Point& x) const {
		Point ret; 
		for (int d=0; d<DIM; d++) {
			ret[d] = (*this)[d] - x[d]; 
		}
		return ret; 
	}

	/// compute distance from a point 
	double Distance(Point x = {0,0,0}) const {
		double dist = 0.; 
		for (int i=0; i<DIM; i++) {
			dist += ((*this)[i] - x[i])*((*this)[i] - x[i]); 
		}
		return sqrt(dist); 
	}

	/// return pointer to data 
	double* GetData() {return _x.data(); }
private:
	/// store location array 
	std::array<double,DIM> _x; 
}; 

std::ostream& operator<<(std::ostream& out, const Point& p); 

} // end namespace fem 