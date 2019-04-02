#pragma once 

#include "General.hpp"
#include "Operator.hpp"
#include "Vector.hpp"

namespace fem 
{

/// base class for all linear solvers 
class Solver {
public:
	/// set the SparseMatrix 
	virtual void SetOperator(const Operator* A) {_A = A; }

	/// solve interface for \f$ A x = b \f$ must be defined 
	/** \param rhs right hand side vector \param x solution vector */ 
	virtual void Solve(Vector& rhs, Vector& x) = 0; 
protected:
	/// left hand side SparseMatrix pointer 
	const Operator* _A; 
	/// number of rows of _A 
	int _m; 
	/// number of columns of _A 
	int _n; 
}; 

} // end namespace fem 