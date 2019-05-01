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
	return 1; 
}

void solve(int nref, int order) {
	HWCounter hwc; 
	hwc.Reset(); 
	SquareMesh mesh(nref, nref, {0,0}, {1,1}); 

	// fem space 
	HWCounter lspace; 
	LagrangeSpace h1(mesh, order); 
	lspace.Read(); 
	lspace.PrintStats("fe space"); 

	// matrix builder 
	HWCounter hwc3; 
	FEMatrix lhs(&h1); 
	lhs.AddIntegrator(new WeakDiffusionIntegrator); 
	lhs.AddIntegrator(new MassIntegrator); 

	// vector builder 
	RHS rhs(&h1); 
	ConstantCoefficient Source(1); 
	rhs.AddIntegrator(new DomainIntegrator(&Source, 0, 1)); 

	// apply the boundary conditions
	lhs.ApplyDirichletBoundary(rhs, 0.);
	lhs.ConvertToBatch();
	hwc3.Read(); 
	hwc3.PrintStats("assembly"); 

	HWCounter hwc2; 
	GridFunction x(&h1); 
	CG cg(&lhs, 1e-10, 1000); 
	hwc2.Reset(); 
	cg.Solve(rhs, x); 
	hwc.Read(); 
	hwc2.Read(); 

	hwc2.PrintStats("cg solve"); 
	hwc.PrintStats(); 
}

int main(int argc, char* argv[]) {
	int N = 20; 
	int p = 1; 
	if (argc > 1) N = atoi(argv[1]); 
	if (argc > 2) p = atoi(argv[2]); 

	solve(N, p); 

	CH_TIMER_REPORT(); 
}