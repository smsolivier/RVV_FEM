#include "FEM.hpp"

using namespace std; 
using namespace fem; 

#define BATCHES 100
#define N 9

bool MVO1(double* mats, int* dofs, Vector& x, Vector& b, Vector& ans) {
	HWCounter hwc; 
#ifdef RV_MVOUTER 
	MVOuter_RV(N, BATCHES, mats, dofs, x.GetData(), b.GetData()); 
#endif
	hwc.Read(); 
	hwc.PrintStats("mvo1"); 

	b -= ans; 
	return b.L2Norm() < 1e-10; 
}

bool MVO2(double* mats, int* dofs, Vector& x, Vector& b, Vector& ans) {
	HWCounter hwc; 
#ifdef RV_MVOUTERC 
	MVOuterC_RV(N, BATCHES, mats, dofs, x.GetData(), b.GetData()); 
#endif
	hwc.Read(); 
	hwc.PrintStats("mvo2"); 

	b -= ans; 
	return b.L2Norm() < 1e-10; 
}

int main() {
#if defined RV_MVOUTER && RV_MVOUTERC 
	double* mats = new double[N*N*BATCHES]; 
	int* dofs = new int[N*BATCHES]; 
	Array<Matrix*> mat_array(BATCHES); 
	Array<Array<int>> dof_array(BATCHES); 

	Vector x(BATCHES*N); 
	Vector b(BATCHES*N); 
	b = 0; 

	for (int i=0; i<x.GetSize(); i++) {
		x[i] = (double)rand()/RAND_MAX;
		// x[i] = 1;  
	}

	for (int n=0; n<BATCHES; n++) {
		dof_array[n].Resize(N); 
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				mats[N*N*n + j + i*N] = n+1; 
			}
			dofs[N*n+i] = N*n + i; 
			dof_array[n][i] = N*n+i; 
		}
		mat_array[n] = new Matrix(N); 
		mat_array[n]->operator=(n+1); 
	}

	Vector ans(BATCHES*N); 
	HWCounter hwc; 
	for (int b=0; b<BATCHES; b++) {
		Vector tmp(N), prod(N); 
		x.GetFromDofs(dof_array[b], tmp); 
		mat_array[b]->Mult(tmp, prod); 
		ans.AddFromDofs(dof_array[b], prod); 
	}
	hwc.Read(); 
	hwc.PrintStats("inner loop vectorization"); 

	TEST(MVO1(mats, dofs, x, b, ans), "outer loop type 1"); 
	b = 0.; 
	TEST(MVO2(mats, dofs, x, b, ans), "outer loop type 2"); 
#endif
}