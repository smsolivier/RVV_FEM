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

void solve(int nref, int order, double& L2, uint64_t& cycles) {
	HWCounter hwc; 
	SquareMesh mesh(nref, nref, {0,0}, {1,1}); 

	// fem space 
	LagrangeSpace h1(mesh, order); 

	// matrix builder 
	FEMatrix lhs(&h1); 
	lhs.AddIntegrator(new WeakDiffusionIntegrator); 
	lhs.AddIntegrator(new MassIntegrator); 

	// vector builder 
	RHS rhs(&h1); 
	FunctionCoefficient Source(source); 
	rhs.AddIntegrator(new DomainIntegrator(&Source, 2, 0)); 

	// apply the boundary conditions
	lhs.ApplyDirichletBoundary(rhs, 0.);
	lhs.ConvertToBatch();

	GridFunction x(&h1); 
	CG cg(&lhs, 1e-10, 1000); 
	cg.Solve(rhs, x); 
	hwc.Read(); 

	FunctionCoefficient ex(exact); 
	L2 = x.L2Error(&ex); 
	cycles = hwc.Cycles(); 
}

int main(int argc, char* argv[]) {
	int N = 8; 
	int p = 1; 
	if (argc > 1) N = atoi(argv[1]); 
	if (argc > 2) p = atoi(argv[2]); 

	Array<int> n = {N, 2*N, 4*N}; 

	double L2; 
	uint64_t cycles; 

	for (int i=0; i<n.GetSize(); i++) {
		solve(n[i], p, L2, cycles); 		
		cout << L2 << " " << cycles << endl; 
	}

}