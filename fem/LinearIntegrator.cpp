#include "LinearIntegrator.hpp"

using namespace std; 

namespace fem 
{

void DomainIntegrator::Assemble(Element& el, Vector& elvec) {
	CH_TIMERS("domain integrator"); 
	Quadrature* quad = QRules.Get(el.GetType(), 
		_oa*el.GetOrder()+_ob, INTEGRATION_TYPE); 
	elvec.SetSize(el.GetNumNodes()); 
	elvec = 0; 
	ElTrans& trans = el.GetTrans(); 
	double c = 1.; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		el.CalcShape(quad->X(n), _shape); 
		if (_c) c = _c->Eval(trans, quad->X(n)); 
		_shape *= quad->Weight(n) * trans.Determinant() * c; 
		elvec += _shape; 
	}
}

void NormalFaceIntegrator::Assemble(Element& el, Vector& elvec) {
	Quadrature* quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	elvec.SetSize(el.GetNumNodes()); 
	elvec = 0.; 
	Vector nor(el.GetMeshDim()); 
	Vector v(el.GetMeshDim()); 
	ElTrans& trans = el.GetTrans(); 
	for (int n=0; n<quad->NumPoints(); n++) {
		Point ip = quad->X(n); 
		trans.SetX(ip); 
		CalcNormal(trans.Jacobian(), nor); 
		_vc->Eval(trans, ip, v); 
		double dot = v * nor; 
		double alpha = .5*(fabs(dot) + dot); 
		double beta = .5*(dot - fabs(dot));  
		el.CalcShape(ip, _shape); 
		double inflow = _inflow->Eval(trans, ip); 
		_shape *= -inflow*beta*quad->Weight(n) * trans.Weight(); 
		elvec += _shape; 
	}
}

void NormalFluxIntegrator::Assemble(Element& el, Vector& elvec) {
	Quadrature* quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	elvec.SetSize(el.GetNumNodes()); 
	elvec = 0.; 
	Vector nor(el.GetMeshDim()); 
	ElTrans& trans = el.GetTrans(); 
	for (int n=0; n<quad->NumPoints(); n++) {
		Point ip = quad->X(n); 
		trans.SetX(ip); 
		CalcNormal(trans.Jacobian(), nor); 
		double c = _c->Eval(ip); 
		el.CalcShape(ip, _shape); 
		_shape *= c * quad->Weight(n); 
		elvec += _shape; 
	}
}

} // end namespace fem 