#include "FEM.hpp"

using namespace std; 
using namespace fem; 

int main(int argc, char* argv[]) {
	int N = 20; 
	int p = 2;
	if (argc>1) N = atoi(argv[1]); 
	if (argc>2) p = atoi(argv[2]);  
	SquareMesh mesh(N, N, {0,0}, {1,1}); 
	LagrangeSpace h1(mesh, p); 
	FEMatrix lhs(&h1); 
	// lhs.AddIntegrator(new WeakDiffusionIntegrator);
	for (int e=0; e<h1.GetNumElements(); e++) {
		lhs[e] = e; 
	}
	GridFunction x(&h1); 
	RHS b(&h1); 
	for (int i=0; i<b.GetSize(); i++) {
		b[i] = (double)rand()/RAND_MAX; 
		x[i] = (double)rand()/RAND_MAX; 
	}	
	lhs.ApplyDirichletBoundary(b, 0.); 
	RHS b2(b); 

	lhs.ConvertToBatch(); 
	HWCounter hwc; 
	lhs.Mult(x, b);
	hwc.Read(); 
	hwc.PrintStats("outer loop matvec");  

	Array<int> vdofs; 
	Vector elvec; 
	Vector prod; 
	int height = (p+1)*(p+1); 
	HWCounter hwc2; 
	for (int e=0; e<h1.GetNumElements(); e++) {
		vdofs = h1.GetVDofs(e);  
		x.GetFromDofs(vdofs, elvec); 
		lhs[e].Mult(elvec, prod); 
		b2.AddFromDofs(vdofs, prod); 
	}
	hwc2.Read(); 
	hwc2.PrintStats("inner loop matvec"); 

	cout << "speedup = " << hwc.FlopsPerCycle()/hwc2.FlopsPerCycle() << endl; 
	cout << "miss ratio = " << hwc.CacheMissRate() / hwc2.CacheMissRate() << endl; 
	cout << "FOM ratio = " << hwc2.GetFOM() / hwc.GetFOM() << endl; 

	b2 -= b; 
	TEST(b2.L2Norm()<1e-10, "fematvec"); 
}