#pragma once 

#include "General.hpp"
#include "Vector.hpp"
#include "Point.hpp"

namespace fem 
{

/// represent a polynomial of one variable 
class Poly1D {
public:
	/// default constructor 
	Poly1D() { }
	/// constructor 
	Poly1D(std::vector<double> c) {SetCoef(c); }
	/// construct with a Vector 
	Poly1D(Vector& c) {
		_c.resize(c.GetSize()); 
		for (int i=0; i<c.GetSize(); i++) {
			_c[i] = c[i]; 
		}
	}

	/// set the coefficients 
	/** evaluates in ascending powers \f[ f(x) = \sum_i c_i x^i \f] */ 
	void SetCoef(std::vector<double> c) {_c = c; }	

	/// evaluate polynomial 
	inline double Eval(double x) const {
		CHECK(_c.size() > 0); 
		double ret = 0; 
		for (int i=0; i<_c.size(); i++) {
			ret += _c[i] * pow(x, i); 
		}
		return ret; 
	}

	/// evaluate polynomial 
	inline double operator()(double x) const {return Eval(x); }

	/// return the derivative of (*this) 
	Poly1D Derivative() const {
		CHECK(_c.size() > 0); 
		std::vector<double> c; 
		for (int i=1; i<_c.size(); i++) {
			c.push_back(_c[i]*i); 
		}
		return Poly1D(c); 
	}

	/// print the coefficients 
	void Print(std::ostream& out=std::cout) const {
		std::string space = ""; 
		if (!EQUAL(_c[0], 0.)) {
			out << _c[0] << " "; 
			space = " "; 
		}
		for (int i=1; i<_c.size(); i++) {
			std::string exp = "^" + std::to_string(i); 
			if (i==1) exp = ""; 
			if (!EQUAL(_c[i], 0.)) {
				std::string sign = "-"; 
				if (_c[i] >= 0) sign = "+"; 
				out << sign << space << fabs(_c[i]) << "x" << exp << " "; 
				space = " "; 
			}
		}
		out << std::endl; 
	}
private:
	/// coefficients of polynomial 
	std::vector<double> _c; 
}; 

/// represent product of 2 or 3 Poly1D's 
class PolyProduct {
public:
	PolyProduct() { }
	/// 1D constructor 
	PolyProduct(const Poly1D& a) {
		_poly.Append(a); 
	}
	/// constructd a 2D polyproduct 
	PolyProduct(const Poly1D& a, const Poly1D& b) {
		_poly.Append(a); 
		_poly.Append(b); 
	}
	/// construct a 3D poly product 
	PolyProduct(const Poly1D& a, const Poly1D& b, const Poly1D& c) {
		_poly.Append(a); 
		_poly.Append(b); 
		_poly.Append(c); 
	}

	/// evaluate 
	inline double operator()(const Point& x) const {
		double ret = 1.;  
		for (int i=0; i<_poly.GetSize(); i++) {
			ret *= _poly[i](x[i]); 
		}
		return ret; 
	}
	inline double Eval(const Point& x) const {return (*this)(x); }

	/// get gradient of this 
	void Gradient(Array<PolyProduct>& grad) const {
		grad.Resize(_poly.GetSize()); 
		if (_poly.GetSize()==1) {
			grad[0] = PolyProduct(_poly[0].Derivative()); 
		} else if (_poly.GetSize()==2) {
			grad[0] = PolyProduct(_poly[0].Derivative(), _poly[1]); 
			grad[1] = PolyProduct(_poly[0], _poly[1].Derivative()); 
		} else if (_poly.GetSize()==3) {
			grad[0] = PolyProduct(_poly[0].Derivative(), _poly[1], _poly[2]); 
			grad[1] = PolyProduct(_poly[0], _poly[1].Derivative(), _poly[2]); 
			grad[2] = PolyProduct(_poly[0], _poly[1], _poly[2].Derivative()); 
		}
	}
private:
	/// store all polynomials 
	Array<Poly1D> _poly; 
}; 

/// multi variate polynomial based on multiplying 1D polynomials 
class Polynomial {
public: 
	/// default constructor 
	Polynomial(); 
	/// constructor 
	Polynomial(std::vector<double> cx, std::vector<double> cy) {
		_c.push_back(Poly1D(cx)); 
		_c.push_back(Poly1D(cy)); 
		_dim = 2; 
	}

	/// evaluate at x in _dim dimensions 
	double Eval(const Point x) const {
		double mult = 1.; 
		for (int i=0; i<_dim; i++) {
			mult *= _c[i].Eval(x[i]); 
		} 
		return mult; 
	}
private:
	/// number of Poly1D copies 
	int _dim; 
	/// coefficients for _dim dimensions 
	std::vector<Poly1D> _c; 
}; 

/// generate Lagrange Polynomials for 1D basis functions 
void GenLagrangePolynomials(int p, double a, double b, Array<Poly1D>& polys); 

} // end namespace fem 