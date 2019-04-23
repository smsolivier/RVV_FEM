#include "CG.hpp" 

using namespace std; 
namespace fem 
{

void CG::Solve(Vector& rhs, Vector& x) {
	CH_TIMERS("cg solve"); 
	CHECK(rhs.GetSize() == _A->Height()); 

	int N = _A->Height(); 
	x.SetSize(N); 

	// set initial guess 
	for (int i=0; i<x.GetSize(); i++) {
		x[i] = 0.; 
	}

	Vector r(N), s(N), As(N), Ax(N); 

	// compute: r = rhs - A*x 
	_A->Mult(x, Ax); 
	Subtract(rhs, Ax, r); 
	// save r 
	s = r; 

	double denom, alpha, beta, norm;

	int iter; 
	chrono::high_resolution_clock::time_point start; 
	chrono::duration<double> elapsed; 
	if (_verbose) cout << "starting CG iterations" << endl; 
	for (iter=0; iter<_max_iter; iter++) {
		if (_verbose) start = chrono::high_resolution_clock::now(); 
		As = 0.; 
		_A->Mult(s, As); 

		// denom = s dot As 
		denom = s.Dot(As); 

		CHECK(denom != 0); 

		// alpha = s dot r/denom 
		alpha = s.Dot(r); 
		alpha /= denom; 

		// x = x + s*alpha 
		Add(x, s*alpha, x); 

		// r = rhs - A*x 
		Ax = 0.; 
		_A->Mult(x, Ax); // Ax = A*x 
		Subtract(rhs, Ax, r); 

		// beta = - r dot As/denom 
		beta = r.Dot(As); 
		beta *= -1./denom; 

		// s = r + s*beta 
		Add(r, s*beta, s); 

		// compute L2 norm of r 
		norm = r.L2Norm(); 

		if (_verbose) {
			elapsed = chrono::high_resolution_clock::now() - start; 
			printf("\titeration %5i, residual = %8.3e, %8.3g s/iter\n", iter, norm, elapsed.count()); 
		}

		if (norm < _tol) break; 
	}

	if (norm < _tol) {
		_converged = true; 
	} else {
		_converged = false; 
		WARNING("maximum number of iterations reached. Final norm = " << norm); 
	}

	if (_print) {
		cout << "number of iterations = " << iter << endl; 
		cout << "final norm = " << norm << endl; 
	}
}

} // end namespace fem 