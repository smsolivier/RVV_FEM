#include "FEM.hpp"

using namespace std; 
using namespace fem; 

int dim = 3; 

double Exact(const Point& x) {
	double sum = 1; 
	for (int i=0; i<dim; i++) {
		sum *= sin(M_PI*x[i]); 
	}
	return sum; 
}

double Source(const Point& x) {
	return 3*M_PI*M_PI*Exact(x); 
}

double ComputeError(int N, int p) {
	HWCounter hwc; 
	CubeMesh mesh({N, N, N}); 
	LagrangeSpace h1(mesh, p); 

	GridFunction x(&h1); 

	FEMatrix K(&h1); 
	K.AddIntegrator(new WeakDiffusionIntegrator); 

	RHS rhs(&h1); 
	FunctionCoefficient source(Source); 
	rhs.AddIntegrator(new DomainIntegrator(&source)); 

	K.ApplyDirichletBoundary(rhs); 

	K.ConvertToBatch(); 
	hwc.Reset(); 
	CG cg(&K, 1e-5, 1000); 
	cg.Solve(rhs, x); 
	hwc.Read(); 

	FunctionCoefficient ex(Exact); 
	double err = x.L2Error(&ex); 

	hwc.PrintStats("N = " + to_string(N)); 
	
	return err; 
}

int main(int argc, char* argv[]) {
	int N = 2; 
	if (argc > 1) N = atoi(argv[1]); 

	double E1, E2, p; 
	E1 = ComputeError(N, 1); 
	E2 = ComputeError(2*N, 1); 

	p = log(E1/E2)/log(2); 
	TEST(abs(p - 2.) < 1e-1, "p = " << p); 
}