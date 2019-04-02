#include "FEM.hpp"

using namespace std; 
using namespace fem; 

int main() {
	int N = 50; 
	Matrix mat(N); 
	mat = 10; 
	bool pass = true; 
	for (int i=0; i<N; i++) {
		if (mat.GetData()[i] != 10) pass = false; 
	}
	TEST(pass, "operator="); 

	pass = true; 
	for (int i=0; i<N; i++) {
		for (int j=0; j<N; j++) {
			mat(i,j) = i*j; 
			if (mat(i,j) != i*j) pass = false; 
		}
	}
	TEST(pass, "index"); 

	Matrix mat2 = mat; 
	pass = true; 
	for (int i=0; i<N*N; i++) {
		if (mat.GetData()[i] != mat2.GetData()[i]) pass = false; 		
	}
	TEST(pass, "copy constructor"); 

	mat.SetSize(N+1, N+1); 
	pass = true; 
	for (int i=0; i<(N+1)*(N+1); i++) {
		if (mat.GetData()[i] != 0.) pass = false; 
	}
	TEST(pass, "resize"); 

	// test inverse 
	Matrix a(2,2); 
	a(0,0) = 2.; a(0,1) = 5.; a(1,1) = .5; 
	Matrix inv; 
	a.Inverse(inv); 
	Matrix eye; 
	a.Mult(inv, eye); 
	TEST(eye(0,0)==1 && eye(1,1)==1 && eye(0,1)==0 && eye(1,0)==0, 
		"2x2 inverse + mult"); 
	Vector X(2); 
	X[0] = 1; X[1] = 2; 
	Vector RHS(2); 
	a.Mult(X, RHS); 
	a.Solve(RHS, X); 
	TEST(EQUAL(X[0], 1.) && EQUAL(X[1], 2.), "solve 2x2"); 

	Matrix m(3,3); 
	m(0,0) = 2.; m(0,1) = 1.; m(1,0) = .5; m(1,1) = 2; m(1,2) = .5; m(2,2) = 1.; 
	m(2,1) = .5; 
	Matrix minv; 
	m.Inverse(minv); 
	Matrix meye; 
	m.Mult(minv, meye); 
	pass = true; 
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			if (i==j) {
				if (abs(meye(i,j) - 1.) > 1e-12) pass = false; 
			} else {
				if (abs(meye(i,j)) > 1e-12) pass = false; 
			}
		}
	}
	TEST(pass, "3x3 inverse + mult"); 

	// test matvec product 
	Vector x(3); 
	for (int i=0; i<3; i++) {
		x[i] = i+1; 
	}
	Vector mx; 
	m.Mult(x, mx); 
	TEST(mx[0]==4 && mx[1]==6 && mx[2]==4, "matvec (3x3)"); 

	// add trans mult 
	Matrix tmult(3); 
	tmult = 1.; 
	m.AddTransMult(meye, tmult); 
	TEST(tmult(0,0)==3 && tmult(0,1)==1.5 && tmult(0,2)==1
		&& tmult(1,0)==2 && tmult(1,1)==3 && tmult(1,2)==1.5
		&& tmult(2,0)==1 && tmult(2,1)==1.5 && tmult(2,2)==2, "add trans mult"); 

	// diagonal check 
	Matrix diag(3,3); 
	diag(0,0) = 1; 
	diag(1,1) = 1; 
	diag(2,2) = 1; 
	TEST(diag.IsDiagonal() && !tmult.IsDiagonal(), "diagonal check is correct"); 

	// blas DGEMM 
	Matrix A(3,3); 
	Matrix B(3,3); 
	Matrix C(3,3); 
	Matrix C2(3,3); 

	for (int i=0; i<A.Height(); i++) {
		for (int j=0; j<A.Width(); j++) {
			A(i,j) = (double)rand()/RAND_MAX; 
			B(i,j) = (double)rand()/RAND_MAX; 
		}
	}

	A.Mult(1., B, 0., C); 

	for (int i=0; i<A.Height(); i++) {
		for (int j=0; j<B.Width(); j++) {
			for (int k=0; k<A.Width(); k++) {
				C2(i,j) += A(i,k)*B(k,j); 
			}
		}
	}
	TEST(C2 == C, "blas dgemm"); 
}