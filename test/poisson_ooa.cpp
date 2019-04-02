#include "FEM.hpp"

using namespace std; 
using namespace fem; 

int dim = 2; 

double exact(const Point& x) {
	double sum = 1; 
	for (int i=0; i<dim; i++) {
		sum *= sin(M_PI*x[i]); 
	}
	return sum; 
}

double source(const Point& x) {
	return (pow(2, dim-1)*M_PI*M_PI + 1)*exact(x); 
}

void solve(int nref, int order, double& LinfE) {
	SquareMesh mesh(nref, nref, {0,0}, {1,1}); 
	// dim = mesh.GetDim(); 

	// fem space 
	LagrangeSpace h1(mesh, order); 

	// matrix builder 
	LHS lhs(&h1); 
	lhs.AddIntegrator(new WeakDiffusionIntegrator); 
	lhs.AddIntegrator(new MassIntegrator); 

	// vector builder 
	RHS rhs(&h1); 
	FunctionCoefficient Source(source); 
	rhs.AddIntegrator(new DomainIntegrator(&Source)); 

	// apply the boundary conditions
	lhs.ApplyDirichletBoundary(rhs, 0.);

	// test symmetry 
	TEST(lhs.IsSymmetric(), "diffusion matrix symmetric"); 

	GridFunction x(&h1); 
	CG cg(&lhs, 1e-10, 1000); 
	cg.Solve(rhs, x); 
	// EigenLU elu(&lhs); 
	// elu.Solve(rhs, x); 
	// EigenBiCGStab solver(&lhs, 1e-10, 1000); 
	// solver.Solve(rhs, x); 

	FunctionCoefficient ex(exact); 
	LinfE = x.L2Error(&ex); 

	// Writer writer; 
	// writer.Add(x, "u"); 
	// writer.Write(); 
}

int main(int argc, char* argv[]) {
	int N = 4; 
	if (argc > 1) N = atoi(argv[1]); 

	double E1, E2; 
	solve(N, 1, E1); 
	solve(2*N, 1, E2); 

	double p = log(E1/E2)/log(2); 
	TEST(abs(p-2.) < 1e-1, "linear FEM p = " << p); 

	N--; 
	solve(N, 2, E1); 
	solve(2*N, 2, E2); 
	p = log(E1/E2)/log(2); 
	TEST(abs(p-3.) < 1e-1, "quadratic FEM p = " << p); 

	N--; 
	solve(N, 3, E1); 
	solve(2*N, 3, E2); 
	p = log(E1/E2)/log(2); 
	TEST(abs(p-4.) < 1e-1, "cubic FEM p = " << p); 
}