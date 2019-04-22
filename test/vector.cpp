#include "FEM.hpp"

using namespace std; 
using namespace fem; 

int main() {
	Vector v(10); 
	v.SetSize(15); 
	TEST(v[14] == 0, "resize"); 

	bool pass = true; 
	for (int i=0; i<v.GetSize(); i++) {
		v[i] = i; 
		if (v[i] != i) pass = false; 
	}
	TEST(pass, "index assignment"); 

	Vector w = v; 
	pass = true; 
	for (int i=0; i<w.GetSize(); i++) {
		if (w[i] != v[i]) pass = false; 
	}
	TEST(pass, "copy constructor"); 
	w[3] = 69; 
	TEST(!EQUAL(w[3], v[3]), "true copy"); 

	Vector a(10); 
	Vector b2(10); 
	for (int i=0; i<10; i++) {
		a[i] = i+1; 
		b2[i] = i; 
	}
	a += b2; 
	for (int i=0; i<10; i++) {
		if (fabs(a[i] - (2*i+1))>1e-10) pass = false; 
	}
	TEST(pass, "vector add"); 

	b2 *= 2.; 
	b2.Print(); 
	for (int i=0; i<10; i++) {
		if (fabs(b2[i] - 2.*i) > 1e-10) pass = false; 
	}
	TEST(pass, "scale vector"); 

#if 0
	// test vector matrix mult 
	Vector x(2); 
	x[0] = 2.; 
	x[1] = 3.; 
	Matrix A(2,2); 
	A(0,0) = 5.; 
	A(0,1) = 9.; 
	A(1,0) = 8.; 
	A(1,1) = 2.; 
	Vector b; 
	x.Mult(A, b); 
	TEST(b[0]==34 && b[1]==24, "vector matrix mult"); 

	Vector y(2); 
	y[0] = 1.; 
	y[1] = 2.; 
	Matrix B(2,3); 
	B(0,0) = 1.; B(0,1) = 2; B(0,2) = 3; 
	B(1,0) = 4; B(1,1) = 5; B(1,2) = 6; 
	y.Mult(B, b); 
	TEST(b[0]==9 && b[1]==12 && b[2]==15, "vector matrix mult (non-square)");

	// test dot product 
	Vector omega(2, 1.); 
	Vector nor(2, -1.); 
	nor[1] = 0.; 
	double dot = nor * omega; 
	TEST(EQUAL(-1., dot), "dot product");
#endif
}