#include "FEM.hpp"

using namespace std; 
using namespace fem; 

#define NRUNS 1
#define N 1000

void VecAdd() {
	Vector a(N); 
	Vector b(N);
	Vector c(N);  

	for (int i=0; i<N; i++) {
		a[i] = 2. ; 
		b[i] = 1.; 
		c[i] = 0; 
	}

	double q = 0; 
	double avl = 0; 
	HWCounter hwc; 
	hwc.Reset();
	for (int n=0; n<NRUNS; n++) {
		Add(a,b,c); 
	}
	hwc.Read();
	hwc.PrintStats("vec add");  
}

void MatVec() {
	Matrix M(N,N); 
	Vector x(N); 
	Vector b(N); 

	HWCounter hwc; 
	M.Mult(x, b); 
	hwc.Read(); 
	hwc.PrintStats("mat vec"); 
}

void BatchMatVec() {
	int n = 4; 
	Array<Matrix*> mats(N); 
	for (int i=0; i<N; i++) {
		mats[i] = new Matrix(n); 
		for (int j=0; j<n*n; j++) {
			mats[i]->operator[](j) = i; 
		}
	}
	Vector x(n); 
	Vector b(n); 

	HWCounter hwc; 
	for (int i=0; i<N; i++) {
		mats[i]->Mult(x, b); 
	}
	hwc.Read(); 
	hwc.PrintStats("batch mat vec"); 
}

void VDot() {
	Vector x(N); 
	Vector y(N); 

	HWCounter hwc; 
	x.Dot(y); 
	hwc.Read(); 
	hwc.PrintStats("vdot"); 
}

int main() {
	VecAdd(); 
	MatVec(); 
	BatchMatVec(); 
	VDot(); 
}