#include "FEM.hpp"

using namespace std; 
using namespace fem; 

int main() {
	HWCounter hwc; 
	int N = 10; 
	Vector a(N); 
	Vector b(N); 

	for (int i=0; i<N; i++) {
		a[i] = 2.; 
		b[i] = 2.; 
	}

	hwc.Reset(); 
	// a += b; 
	a *= 2.; 
	hwc.Read(); 

	a.Print(); 
	// cout << hwc.AvgVecLen() << endl; 
}