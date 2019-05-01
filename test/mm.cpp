#include "FEM.hpp"

using namespace std; 
using namespace fem; 

bool MatMult() {
	int N = 14; 
	int M = 18; 
	int K = 7; 
	Matrix A(N,M); 
	Matrix B(M,K); 
	Matrix C(N,K); 
	Matrix ans(N,K); 

	for (int i=0; i<N; i++) {
		for (int j=0; j<M; j++) {
			A(i,j) = (double)rand()/RAND_MAX; 
			// A(i,j) = i; 
		}
	}

	for (int i=0; i<M; i++) {
		for (int j=0; j<K; j++) {
			B(i,j) = (double)rand()/RAND_MAX; 
			// B(i,j) = 2; 
		}
	}

	for (int i=0; i<N; i++) {
		for (int j=0; j<K; j++) {
			C(i,j) = 0; 
			ans(i,j) = 0; 
		}
	}

	A.Mult(B,C); 

	for (int i=0; i<N; i++) {
		for (int j=0; j<K; j++) {
			for (int k=0; k<M; k++) {
				ans(i,j) += A(i,k) * B(k,j); 
			}
		}
	}
	C -= ans; 
	return C.FrobeniusNorm() < 1e-10; 
}

bool ATM() {
	int N = 2; 
	int M = 4; 
	Matrix A(N,M); 
	Matrix B(N,M); 
	Matrix C(M,M); 

	for (int i=0; i<N; i++) {
		for (int j=0; j<M; j++) {
			A(i,j) = 1; 
			B(i,j) = 2; 
		}
	}

	C = -1.; 
	A.AddTransMult(B, C); 
	C.Print(); 

	return true; 
}

int main() {
	TEST(MatMult(), "matmult"); 
	// TEST(ATM(), "Add Trans Mult"); 
	ATM(); 
}