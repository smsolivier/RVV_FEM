#pragma once 

#include "General.hpp"
#include "Quadrature.hpp"
#include "ElTrans.hpp"
#include "Coefficient.hpp"

namespace fem 
{

/// represent an integrator for the LHS 
class BilinearIntegrator {
public:
	/// constructor 
	BilinearIntegrator() { }
	/// assemble the local matrix 
	virtual void Assemble(Element& el, Matrix& elmat) {
		ERROR("not implemented"); 
	}
	/// assemble the local matrix for a mixed system 
	/** \param[in] trial trial space element (solution is expanded in trial space) 
		\param[in] test test space element (ie space used for weak formulation) 
		\param[out] elmat assembled local matrix 
	*/ 
	virtual void MixedAssemble(Element& trial, Element& test, 
		Matrix& elmat) {
		ERROR("not implemented"); 
	}
	/// assemble a face integrator 
	/** \param[in] e current element 
		\param[in] ep neighboring element 
		\param[in] face element transformation describing the face between e and ep 
		\param[out] elmat assembled local face matrix 
	*/ 
	virtual void AssembleFaceMatrix(Element* e, Element* ep, 
		FaceTransformations& fts, Matrix& elmat) {
		ERROR("not implemented"); 
	}
protected:
};

/// integrate \f$ \int \nabla B_i \cdot \nabla B_j dV \f$ 
class WeakDiffusionIntegrator : public BilinearIntegrator {
public:
	/// constructor 
	WeakDiffusionIntegrator(Coefficient* c=NULL) {_c = c; }
	/// assemble 
	void Assemble(Element& el, Matrix& elmat); 
private:
	/// store grad shape matrix in physical space 
	Matrix _pgshape; 
	/// temporary matrix 
	Matrix _tmp; 
	/// coefficient 
	Coefficient* _c; 
}; 

/// integrate \f$ \int B_i B_j dV \f$ 
class MassIntegrator : public BilinearIntegrator {
public:
	/// constructor 
	MassIntegrator(Coefficient* c=NULL) {_c = c; }
	/// assemble 
	void Assemble(Element& el, Matrix& elmat); 
	/// assemble mixed system 
	void MixedAssemble(Element& trial, Element& test, 
		Matrix& elmat); 
private:
	/// store shape functions 
	Vector _shape; 
	/// store shape function evaluation in mixed case 
	Vector _shape2; 
	/// store outer product 
	Matrix _op; 
	/// coefficient 
	Coefficient* _c; 
}; 

/// integrate \f$ \int B_i B_j dV \f$ with trapezoidal quadrature 
class MassLumpingIntegrator : public BilinearIntegrator {
public:
	/// constructor 
	MassLumpingIntegrator(Coefficient* c=NULL) {_c = c; }
	/// assemble 
	void Assemble(Element& el, Matrix& elmat); 
private:
	/// store shape functions 
	Vector _shape; 
	/// store outer product 
	Matrix _op; 
	/// store the coefficient 
	Coefficient* _c; 
}; 

/// integrate \f$ \int \vec{S}_i \cdot \vec{S}_j dV \f$ 
class VectorMassIntegrator : public BilinearIntegrator {
public:
	VectorMassIntegrator(Coefficient* c=NULL) {_c = c; }
	/// assemble 
	void Assemble(Element& el, Matrix& elmat); 
private:
	/// store shape functions 
	Vector _shape; 
	/// store outer product 
	Matrix _op; 
	/// constant multiplier 
	Coefficient* _c; 
}; 

/// integrate \f$ \int B_i \vec{v} \cdot \nabla B_j \f$ 
class ConvectionIntegrator : public BilinearIntegrator {
public:
	ConvectionIntegrator(const VectorCoefficient* vc) {_vc = vc; }
	/// assemble 
	void Assemble(Element& el, Matrix& elmat); 
private:
	/// store coefficient 
	const VectorCoefficient* _vc; 
	/// store grad shape matrix 
	Matrix _pgshape; 
	/// store shape evaluation 
	Vector _shape; 
	/// store function evaluation 
	Vector _v; 
	/// intermediate matrix 
	Vector _tmp; 
	/// store outer product 
	Matrix _op; 
}; 

/// integrate \f$ \int B_i \nabla \cdot \vec{S}_j dV \f$ 
class VectorDivergenceIntegrator : public BilinearIntegrator {
public:
	/// constructor 
	VectorDivergenceIntegrator() { }
	/// assemble 
	// void Assemble(Element& el, Matrix& elmat); 
	/// assemble mixed system 
	void MixedAssemble(Element& trial, Element& test, Matrix& elmat); 
private:
	/// store trial shape eval 
	Vector _shape; 
	/// store test grad shape eval 
	Matrix _gshape; 
	/// store divshape 
	Vector _divshape; 
	/// store the outer product 
	Matrix _op; 
}; 

/// integrate face with upwinding 
/** integrates \f$ \int_{\partial V} \frac{1}{2}(| \hat{\Omega} \cdot \hat{n} | 
	+ \hat{\Omega} \cdot \hat{n}) \psi_e + 
	\frac{1}{2}(| \hat{\Omega} \cdot \hat{n} | + 
	\hat{\Omega} \cdot \hat{n}) \psi_{e'} \f$ 
	where \f$ \psi_{e} \f$ is the current element and 
	\f$ \psi_{e'} \f$ is the upwind neighbor 
*/ 
class UpwindFaceIntegrator : public BilinearIntegrator {
public:
	/// constructor 
	UpwindFaceIntegrator(const VectorCoefficient* vc) {_vc = vc; }
	/// assemble face matrix 
	void AssembleFaceMatrix(Element* e, Element* ep, 
		FaceTransformations& fts, Matrix& elmat); 
private:
	/// velocity vector coefficient 
	const VectorCoefficient* _vc; 
	/// store shape functions of element e
	Vector _shape; 
	/// shape functions for element ep 
	Vector _pshape; 
	/// store downwind sub matrix 
	Matrix _d; 
	/// upwind sub matrix 
	Matrix _u; 
}; 

} // end namespace fem 