#pragma once 

#include "General.hpp"
#include "Vector.hpp"
#include "FESpace.hpp"
#include "Quadrature.hpp"
#include "Coefficient.hpp"

namespace fem 
{

/// represent a solution variable on an FESpace 
class GridFunction : public Vector {
public:
	/// default constructor 
	GridFunction(); 
	/// constructor 
	GridFunction(const FESpace* space); 
	/// set the FESpace 
	void SetSpace(const FESpace* space); 
	/// return the FESpace pointer 
	const FESpace* GetSpace() const {return _space; }
	/// evaluate a function at the FEM nodes 
	void Project(double (*f)(const Point&)); 
	/// return the energy norm 
	double EnergyNorm(Quadrature* quad=NULL) const; 
	/// compute L2 error from a coefficient (through integration) 
	double L2Error(Coefficient* c) const; 

	using Vector::operator=; 
protected:
	/// pointer to FESpace object 
	const FESpace* _space; 
}; 

} // end namespace fem 