#pragma once 

#include "General.hpp"
#include "Operator.hpp"
#include "Vector.hpp"
#include "FESpace.hpp"
#include "Array.hpp"
#include "Quadrature.hpp"
#include "ElTrans.hpp"
#include "BilinearIntegrator.hpp"
#include "RHS.hpp"
#include "SparseMatrix.hpp"

#ifdef USE_RISCV
// b batches of NxN scattered matvecs 
extern "C" void MVOuter_RV(int N, int B, const double* mats, 
	const int* dofs, const double* x, double* b); 
// b batches of 4x4 matrices 
extern "C" void MVOuter4_RV(int N, int B, const double* mats, 
	const int* dofs, const double* x, double* b); 
// b batches of 9x9 matrices 
extern "C" void MVOuter9_RV(int N, int B, const double* mats, 
	const int* dofs, const double* x, double* b); 
// b batches of 16x16 matrices 
extern "C" void MVOuter16_RV(int N, int B, const double* mats, 
	const int* dofs, const double* x, double* b);
// b batches of matvecs outputting to contiguous b 
extern "C" void MVOuterC_RV(int N, int B, const double* mats, 
	const int* dofs, const double* x, double* b); 
// scatter add with reduced vl to avoid aliasing 
extern "C" void BatchAdd_RV(int N, int B, const int* dofs, double* ball, double* b); 
#endif

namespace fem 
{

/// store an element by element finite element matrix 
class FEMatrix : public Operator {
public:
	/// default constructor 
	FEMatrix() { }
	/// deconstructor 
	~FEMatrix(); 
	/// construct and set size 
	FEMatrix(const FESpace* space); 
	/// matrix vector product 
	void Mult(const Vector& x, Vector& b) const; 
	/// access elemental matrices 
	Matrix& operator[](int el) {return *_data[el]; }
	/// const access to elemental matrices 
	const Matrix& operator[](int el) const {
		CH_TIMERS("get el matrix"); 
		return *_data[el]; 
	}
	/// add a bilinear integrator 
	void AddIntegrator(BilinearIntegrator* integ); 
	/// apply dirichlet boundary conditions 
	void ApplyDirichletBoundary(RHS& rhs, double val=0); 
	/// convert from element by element to a general SparseMatrix 
	void ConvertToSparseMatrix(SparseMatrix& spmat) const; 
	/// extract the diagonal of the assembled matrix 
	void GetDiagonal(Vector& diag) const; 
	/// diagonal preconditioning on this matrix 
	void DiagonalPrecondition(const Vector& diag); 
	/// convert to batch storage format 
	void ConvertToBatch(); 

	/// subtract two FEMatrix's 
	void operator-=(const FEMatrix& A); 
protected:
	/// store the elemental matrices 
	Array<Matrix*> _data; 
	/// store the fespace 
	const FESpace* _space; 
	/// store data for all matrices contiguously 
	Array<double> _mats; 
	/// store the vdofs for each element contiguously 
	Array<int> _vdofs; 
}; 

} // end namespace fem 