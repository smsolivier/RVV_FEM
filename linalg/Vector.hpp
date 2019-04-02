#pragma once 

#include "General.hpp"
#include "Matrix.hpp"
#include "Array.hpp"

namespace fem 
{

class Matrix; 

/// represent an array 
class Vector : public Array<double> {
public:
	/// default constructor 
	Vector(); 
	/// constructor 
	/** \param len initial size of array 
		\param val initialized value 
	*/ 
	Vector(int len, double val=0); 
	/// get a Vector of values from this corresponding to vdofs 
	void GetFromDofs(const Array<int>& vdofs, Vector& v) const; 
	/// add to this from corresponding vdofs 
	void AddFromDofs(const Array<int>& vdofs, Vector& v); 
	/// resize the array 
	void SetSize(int size, double val=0) {Resize(size); } 
	/// scale by a scalar constant f
	Vector operator*(double val) const; 
	/// scale this by val 
	void operator*=(double val); 
	/// add a Vector to *this 
	void operator+=(const Vector& a); 
	/// subtract a Vector from *this 
	void operator-=(const Vector& a); 
	/// divide this by a vector 
	void operator/=(const Vector& a); 
	/// multiply this by a vector 
	void operator*=(const Vector& a); 
	/// outer product 
	void OuterProduct(const Vector& a, Matrix& b) const; 
	/// vector matrix multiplication 
	void Mult(const Matrix& a, Vector& b) const; 
	/// flip array 
	void Transpose(); 
	/// take the square root of all entries 
	void SquareRoot();

	/// dot product of two vectors 
	double Dot(const Vector& x) const; 
	/// dot product of two vectors 
	double operator*(const Vector& x) const {return Dot(x); }
	/// L2 norm of this 
	double L2Norm() const; 
	/// L infinity norm of this 
	double LinfNorm() const; 

	/// generate N equally spaced points between a and b 
	void Linspace(int N, double a=0, double b=1); 

	/// print to stream 
	void Print(std::ostream& out=std::cout, std::string end="\n") const; 

	using Array<double>::operator=; 
private:
};

/// add two vectors \f$ c = a+b \f$ 
void Add(const Vector& a, const Vector& b, Vector& c);  
/// subtract two vectors \f$ c = a - b \f$ 
void Subtract(const Vector& a, const Vector& b, Vector& c); 
/// generalized add \f$ alpha*a + beta*b = c \f$ 
void Add(double alpha, const Vector& a, double beta, const Vector& b, Vector& c); 

} // end namespace fem 