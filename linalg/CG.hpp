#pragma once 

#include "General.hpp" 
#include "SparseMatrix.hpp" 
#include "IterativeSolver.hpp"

namespace fem 
{

/// Conjugate Gradient solver 
class CG : public IterativeSolver {
public: 
	/// default constructor 
	CG() : IterativeSolver() { }
	/// constructor 
	CG(const Operator* A, double tol=1e-6, int max_iter=25, bool verbose=false)
		: IterativeSolver(A, tol, max_iter, verbose) {	}

	/// solve interface 
	void Solve(Vector& rhs, Vector& x); 
};

} // end namespace fem 
