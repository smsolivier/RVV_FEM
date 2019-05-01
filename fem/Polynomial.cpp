#include "Polynomial.hpp"

namespace fem 
{

void GenLagrangePolynomials(int p, double a, double b, Array<Poly1D>& polys) {
	CH_TIMERS("generate lagrange polynomials"); 
	int N = p + 1; 
	Vector x; 
	x.Linspace(N, a, b); 

	Array<Vector> coefs; 

	for (int k=0; k<N; k++) {
		Matrix A(N); 
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				A(i,j) = pow(x[i], j); 
			}
		}
		Vector b(N); 
		b[k] = 1.; 
		Vector c; 
		A.Solve(b, c); 
		coefs.Append(c); 
	}

	polys.Resize(N); 
	for (int i=0; i<N; i++) {
		polys[i] = Poly1D(coefs[i]); 
	}
}

} // end namespace fem 