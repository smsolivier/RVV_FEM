#pragma once 

#include "General.hpp" 
#include "Solver.hpp"

namespace fem 
{

/// abstract class for iterative linear solvers 
class IterativeSolver : public Solver {
public:
	/// default constructor 
	IterativeSolver() {
		_converged = false; 
		_verbose = false; 
		_print = false; 
		_A = NULL; 
		_tol = 1e-6; 
		_max_iter = 50; 
	}
	/// constructor. store a reference to the SparseMatrix 
	IterativeSolver(const Operator* A, double tol=1e-6, int max_iter=50, 
		bool verbose=false) {
		_A = A; 
		_tol = tol; 
		_max_iter = max_iter; 
		_converged = false; 
		_verbose = verbose; 
		_print = false; 
	}
	/// set the matrix 
	void SetOperator(const Operator* A) {_A = A; }
	/// set the tolerance 
	void SetTol(double tol) {_tol = tol; }
	/// set maximum number of iterations 
	void SetMaxIter(int max_iter) {_max_iter = max_iter; }

	/// solve interface from Solver 
	virtual void Solve(Vector& rhs, Vector& x) = 0; 

	/// return if previous call to solve converged 
	bool GetConverged() const {return _converged; } 
	/// set verbosity 
	void SetVerbose() {_verbose = true; }
	/// print stats after solve 
	void PrintStats() {_print = true; }
protected:
	/// true if previous Solve call converged 
	bool _converged; 
	/// solver tolerance 
	double _tol; 
	/// maximum number of iterations 
	int _max_iter; 
	/// print statistics if true 
	bool _print; 
	/// print residual at each iteration 
	bool _verbose; 
}; 

} // end namespace fem 