#pragma once 

#include "General.hpp"
#include "Vector.hpp"
#ifdef USE_EIGEN 
#include "Sparse"
#endif

namespace fem 
{

/// abstract class for operators 
class Operator {
public:
	/// default constructor set sizes to zero 
	Operator() {
		_m = 0; 
		_n = 0; 
	}
	/// constructor
	/** \param M number of rows 
		\param N number of columns 
	*/ 
	Operator(int M, int N=-1) {Resize(M, N); }

	/// reset the size of the operator 
	void Resize(int M, int N=-1) {
		_m = M; 
		if (N>0) {
			_n = N; 
		} else {
			_n = _m; 
		}
	}

	/// application of operator onto a vector 
	virtual void Mult(const Vector& x, Vector& b) const = 0; 

	#ifdef USE_EIGEN 
	virtual void GetEigenFormat(Eigen::SparseMatrix<double>& eigen) const {
		ERROR("conversion to Eigen not supported"); 
	}
	#endif
	/// get the number of rows 
	int Height() const {return _m; }
	/// get the number of columns 
	int Width() const {return _n; }
protected:
	/// number of rows 
	int _m; 
	/// number of columns 
	int _n; 
}; 

} // end namespace fem 