#pragma once 

#include "General.hpp"
#include "Vector.hpp"
#include "Array.hpp"

#ifdef USE_RISCV
// matrix vector product 
extern "C" void MatVec_RV(int N, int M, const double* A, const double* x, double* b); 
// dgemm: C = A*B 
extern "C" void MatMult_RV(int N, int M, const double* A, const double* B, double* C); 
// scalar mat vec 
extern "C" void MatVec_S(int N, int M, const double* A, const double* x, double* b); 
#endif

namespace fem 
{

class Vector; 

/// represent a dense matrix 
class Matrix : public Array<double> {
public:
	/// default constructor 
	Matrix(); 
	/// constructor (defaults to square matrix) 
	/** \param m number of rows \param n number of columns */ 
	Matrix(int m, int n=-1); 
	/// copy constructor 
	Matrix(const Matrix& m); 
	/// copy assignment 
	Matrix& operator=(const Matrix& m); 

	/// access to data 
	double& operator()(int i, int j); 
	/// const access to data 
	double operator()(int i, int j) const; 
	/// resize the matrix 
	void SetSize(int m, int n=-1); 
	/// return the number of rows 
	int Height() const {return _m; }
	/// return the number of columns 
	int Width() const {return _n; }
	/// set all values to val 
	// void operator=(double val); 
	/// scale all entries by val 
	void operator*=(double val); 
	/// add a submatrix to this 
	void AddMatrix(double a, const Matrix& mat, int i, int j); 
	/// AddMatrix with \f$ a = 1 \f$ 
	void AddMatrix(const Matrix& mat, int i, int j) {AddMatrix(1., mat, i, j); }
	/// add a matrix and store in this 
	void operator+=(const Matrix& a); 
	/// subtract a matrix from this 
	void operator-=(const Matrix& a); 
	/// compute inverse 
	void Inverse(Matrix& inv) const; 
	/// return determinant 
	double Determinant() const; 
	/// return the "weight" (basically the determinant but generalized non-square) 
	double Weight() const; 
	/// matrix matrix multiplication \f$ b = (*this)*a \f$ 
	void Mult(const Matrix& a, Matrix& b) const; 
	/// BLAS matrix matrix product \f$ C = \alpha A B + \beta C \f$ 
	void Mult(double alpha, const Matrix& B, double beta, Matrix& C) const; 
	/// transpose(matrix) matrix multiplication \f$ b += (*this)^T * a \f$ 
	void AddTransMult(const Matrix& a, Matrix& b) const; 
	/// matrix vector product (does not accumulate into b) 
	void Mult(const Vector& x, Vector& b) const; 
	/// Lapack matrix vector product \f$ b = alpha this x + beta * b \f$ 
	void Mult(double alpha, const Vector& x, double beta, Vector& b) const; 
	/// check if *this is symmetric 
	bool IsSymmetric() const; 
	/// check if two matrices are equal 
	bool operator==(const Matrix& mat) const; 
	/// convert gradients of shape functions to divergence vector 
	void GradToDiv(Vector& divshape) const; 
	/// solve a dense system with Lapack 
	void Solve(const Vector& b, Vector& x) const; 
	/// tranpose the matrix 
	void Transpose(Matrix& trans) const; 

	/// compute frobenius norm (square root of sum of squares of all entries) 
	double FrobeniusNorm() const; 

	/// check if matrix is diagonal 
	bool IsDiagonal() const; 

	/// print matrix 
	void Print(std::ostream& out=std::cout) const; 

	using Array<double>::operator=; 
private:
	/// invert a 2x2 matrix directly 
	void Inverse_2x2(Matrix& inv) const; 
	/// invert a 3x3 matrix directly 
	void Inverse_3x3(Matrix& inv) const; 
	/// gaussian elimination 
	int GaussElim(int dim, Vector& x, const Vector& b) const; 
	/// number of rows 
	int _m; 
	/// number of columns 
	int _n; 
}; 

void CalcNormal(const Matrix& J, Vector& nor, bool normalize=true); 

} // end namespace fem 