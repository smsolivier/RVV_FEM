#pragma once 

#include "General.hpp"
#include "Quadrature.hpp"
#include "ElTrans.hpp"
#include "FESpace.hpp"
#include "BilinearIntegrator.hpp"
#include "SparseMatrix.hpp"
#include "RHS.hpp"

namespace fem 
{

/// build the left hand side matrix 
class LHS : public SparseMatrix {
public:
	/// constructor 
	LHS(const FESpace* space, Quadrature* gq=NULL); 
	/// add an integrator 
	void AddIntegrator(BilinearIntegrator* integ); 
	/// assemble a local face integrator 
	void AddFaceIntegrator(BilinearIntegrator* integ); 
	/// apply dirichlet boundary conditions by 
	/// eliminating rows corresponding to boundary nodes 
	void ApplyDirichletBoundary(RHS& rhs, double val=0); 
	/// return the FESpace pointer 
	const FESpace* GetSpace() const {return _space; }
private:
	/// FESpace 
	const FESpace* _space; 
	/// number of elements in FESpace 
	int _nel; 
	/// number of nodes in FESpace 
	int _nnodes; 
	/// quadrature class 
	Quadrature* _gq; 
}; 

/// build the left hand side sparse matrix for a mixed FEM system 
class MixedLHS : public SparseMatrix {
public:
	/// constructor 
	MixedLHS(const FESpace* trial, const FESpace* test, Quadrature* quad=NULL); 
	/// add an integrator 
	void AddIntegrator(BilinearIntegrator* integ); 

	/// return the trial FESpace 
	const FESpace* GetTrialSpace() const {return _trial; }
	/// return the test FESpace 
	const FESpace* GetTestSpace() const {return _test; }
private:
	/// trial space 
	const FESpace* _trial; 
	/// test space 
	const FESpace* _test; 
	/// quadrature class 
	Quadrature* _quad; 
	/// number of elements in spaces 
	int _nel; 
}; 

} // end namespace fem 