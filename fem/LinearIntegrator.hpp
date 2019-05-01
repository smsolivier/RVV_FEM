#pragma once 

#include "General.hpp"
#include "Quadrature.hpp"
#include "ElTrans.hpp"
#include "Coefficient.hpp"
#include "Vector.hpp"

namespace fem 
{

/// abstract class for RHS construction
class LinearIntegrator {
public:
	/// constructor 
	LinearIntegrator() { }
	/// assemble the contribution from el 
	virtual void Assemble(Element& el, Vector& elvec) {
		ERROR("not implemented"); 
	}
protected:
}; 

/// integrate \f$ \int B_i Q dV \f$ 
class DomainIntegrator : public LinearIntegrator {
public:
	/// constructor 
	DomainIntegrator(Coefficient* c=NULL, int oa=2, int ob=1) {
		_c = c; 
		_oa = oa; 
		_ob = ob; 
	}
	/// assemble local vector 
	void Assemble(Element& el, Vector& elvec); 
private:
	/// store shape functions 
	Vector _shape; 
	/// store function for Q 
	Coefficient* _c; 
	/// quadrature order function coefficient 
	int _oa; 
	int _ob; 
}; 

/// integrate boundary contribution for UpwindFaceIntegrator 
class NormalFaceIntegrator : public LinearIntegrator {
public:
	/// constructor 
	NormalFaceIntegrator(VectorCoefficient* vc, Coefficient* inflow) {
		_vc = vc; 
		_inflow = inflow; 
	}
	/// assemble a local vector 
	void Assemble(Element& el, Vector& elvec); 
private:
	/// vector coefficient for angle 
	VectorCoefficient* _vc; 
	/// coefficient for the inflow conditions 
	Coefficient* _inflow; 
	/// store shape evaluations 
	Vector _shape; 
	/// store outer product 
	Matrix _op; 
}; 

/// normal flux integrator 
class NormalFluxIntegrator : public LinearIntegrator {
public:
	/// constructor 
	NormalFluxIntegrator(Coefficient* c) {_c = c; }
	/// assemble a local vector 
	void Assemble(Element& el, Vector& elvec); 
private:
	/// coefficient for the normal flux 
	Coefficient* _c; 
	/// store the shape evaluations 
	Vector _shape; 
}; 

} // end namespace fem 