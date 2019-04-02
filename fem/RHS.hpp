#pragma once 

#include "General.hpp"
#include "Quadrature.hpp"
#include "ElTrans.hpp"
#include "FESpace.hpp"
#include "Vector.hpp"
#include "LinearIntegrator.hpp"

namespace fem 
{

/// build the right hand side vector 
class RHS : public Vector {
public:
	/// constructor 
	RHS(const FESpace* space, Quadrature* gq=NULL); 
	/// add an integrator 
	void AddIntegrator(LinearIntegrator* integ); 
	/// add boundary integrator 
	void AddBoundaryIntegrator(LinearIntegrator* integ); 

	using Vector::operator=; 
protected:
	/// pointer to FESpace 
	const FESpace* _space; 
	/// number of nodes 
	int _nnodes; 
}; 

} // end namespace fem 