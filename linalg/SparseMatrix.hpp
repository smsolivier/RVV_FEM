#pragma once 

#include "General.hpp"
#include "Operator.hpp"
#include "Vector.hpp"
#ifdef USE_EIGEN
#include "Sparse"
#endif

namespace fem 
{

/// store a matrix in compressed column format 
class SparseMatrix : public Operator {
public:
	/// default constructor 
	SparseMatrix(); 
	/// construct and set size 
	SparseMatrix(int M, int N=-1); 
	/// resize the extents of the matrix 
	void Resize(int M, int N=-1); 
	/// copy constructor 
	SparseMatrix(const SparseMatrix& sp); 

	/// matrix vector product \f$ \mathbf{A} x += b \f$ 
	void Mult(const Vector& x, Vector& b) const; 
	/// scale all elements by val 
	void operator*=(double val); 
	
	/// access elements of the matrix 
	/** returns  a a reference to the entry \n
		if the entry does not exist one is created and set to zero \n
		only use if adding non-zeros otherwise it will clutter the sparsematrix 
	*/ 
	double& operator()(int row, int col); 
	/// same as operator() 
	double& Get(int row, int col) {return (*this)(row, col); }

	/// const access to the matrix 
	/** return zeros if entry does not exist but doesn't create new entries */ 
	double operator()(int row, int col) const; 
	/// same as const operator() 
	double At(int row, int col) const; 
	/// remove values less than a given tolerance 
	void ClearNonZeros(double tol=1e-12); 
	/// subtract column rc from rhs and set row rc to zero except for one on the diagonal
	void EliminateRowIntoRHS(int rc, Vector& rhs, double val); 

	/// print a dense matrix 
	void Print(std::ostream& stream = std::cout) const; 

	/// output sparsity information 
	void Sparsity(std::ostream& stream = std::cout) const; 

	/// return number of non-zeros 
	int GetNNZ() const {return _nnz; }
	/// tests if *this is symmetric 
	bool IsSymmetric() const; 
	/// get transpose of matrix 
	void Transpose(SparseMatrix& transpose) const; 

	/// return val, row, colptr format for superlu format 
	void GetSLUFormat(double* val, int* row, int* colptr) const; 
	/// get Eigen matrix 
#ifdef USE_EIGEN 
	void GetEigenFormat(Eigen::SparseMatrix<double>& eigen) const; 
#endif
private:
	/// store a zero to give as a reference 
	double _zero; 
	/// number of non-zero entries 
	int _nnz; 
	/// store the non-zero data 
	std::vector<std::vector<double>> _data; 
	/// store where the non-zeros are 
	std::vector<std::vector<int>> _rowIndex;  
}; 

} // end namespace fem 