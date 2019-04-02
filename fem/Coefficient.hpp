#pragma once 

#include "General.hpp"
#include "ElTrans.hpp"

namespace fem 
{

/// abstract class for evaluating functions
class Coefficient {
public:
	/// evaluate the coefficient at a point in physical space 
	virtual double Eval(const Point& phys_x) const {
		ERROR("base class is not callable. Use pointer to call derived class's Eval"); 
	}
	/// evaluate by first transforming to physical space 
	double Eval(ElTrans& trans, const Point& ref_x) const {
		Point phys; 
		trans.Transform(ref_x, phys); 
		return Eval(phys); 
	}
}; 

/// store a constant value 
class ConstantCoefficient : public Coefficient {
public:
	/// constructor 
	ConstantCoefficient(double c) {_c = c; }
	/// evaluate 
	double Eval(const Point& phys_x) const {return _c; }
private:
	/// constant value 
	double _c; 
}; 

/// store a function of space 
class FunctionCoefficient : public Coefficient {
public:
	/// constructor 
	FunctionCoefficient(double (*f)(const Point& x)) {_f = f; }
	/// evaluate 
	double Eval(const Point& x) const {
		return _f(x); 
	}
private:
	/// store the function pointer 
	double (*_f)(const Point&); 
}; 

/// abstract class for evaluating functions that return vectors 
class VectorCoefficient {
public:
	/// evaluate the vector function in physical space 
	/** \param[in] x_phys point in physical space 
		\param[out] vector output 
	*/ 
	virtual void Eval(const Point& x_phys, Vector& v) const {
		ERROR("not implemented"); 
	}
	/// evaluate in reference space 
	void Eval(ElTrans& trans, const Point& x_ref, Vector& v) const {
		Point x_phys; 
		trans.Transform(x_ref, x_phys); 
		Eval(x_phys, v); 
	}
}; 

/// constant vector coefficient 
class ConstantVectorCoefficient : public VectorCoefficient {
public:
	/// constructor 
	ConstantVectorCoefficient(const Vector& v) {_v = v; }
	/// evaluate 
	void Eval(const Point& x_phys, Vector& v) const {
		v = _v; 
	}
private:
	/// store vector 
	Vector _v; 
}; 

/// evaluate vector functions 
class VectorFunctionCoefficient : public VectorCoefficient {
public:
	/// constructor 
	VectorFunctionCoefficient(void (*f)(const Point&, Vector&)) {_f = f; }
	/// evaluate 
	void Eval(const Point& x, Vector& v) const {
		_f(x, v); 
	}
private:
	/// store vector function 
	void (*_f)(const Point&, Vector&); 
}; 

/// evaluate the product of two coefficients as one coefficient 
class ProductCoefficient : public Coefficient {
public:
	/// constructor 
	ProductCoefficient(Coefficient* c1, Coefficient* c2) {_c1 = c1; _c2 = c2; }
	/// evaluate c1 * c2 
	double Eval(const Point& x) const {
		return _c1->Eval(x) * _c2->Eval(x); 
	}
private:
	/// store the coefficients 
	Coefficient* _c1; 
	/// store the coefficients 
	Coefficient* _c2; 
}; 

} // end namespace fem 