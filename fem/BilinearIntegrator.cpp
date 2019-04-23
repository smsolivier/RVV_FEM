#include "BilinearIntegrator.hpp"

using namespace std; 

namespace fem 
{

void WeakDiffusionIntegrator::Assemble(Element& el, Matrix& elmat) {
	CH_TIMERS("weak diffusion assemble"); 
	Quadrature* quad = QRules.Get(el.GetType(), el.GetOrder()+1, INTEGRATION_TYPE); 
	ElTrans& trans = el.GetTrans();
	elmat.SetSize(el.GetNumNodes()); 
	elmat = 0; 
	_pgshape = 0; 
	_tmp = 0; 
	double c = 1.; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		el.CalcPhysGradShape(trans, _pgshape); 
		_tmp = _pgshape; 
		if (_c) c = _c->Eval(trans, quad->X(n)); 
		_tmp *= quad->Weight(n) * trans.Determinant() * c; 
		_pgshape.AddTransMult(_tmp, elmat); 
	}
	CHECK(elmat.IsSymmetric()); 
}

void MassIntegrator::Assemble(Element& el, Matrix& elmat) {
	CH_TIMERS("mass assemble"); 
	Quadrature* quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	ElTrans& trans = el.GetTrans();
	elmat.SetSize(el.GetNumNodes()); 
	elmat = 0; 
	double c = 1.; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		el.CalcShape(quad->X(n), _shape);
		_shape.OuterProduct(_shape, _op);
		if (_c) c = _c->Eval(trans, quad->X(n)); 
		_op *= c * quad->Weight(n) * trans.Determinant(); 
		elmat += _op;
	}
}

void MassIntegrator::MixedAssemble(Element& trial_fe, Element& test_fe, 
	Matrix& elmat) {
	Quadrature* quad = QRules.Get(trial_fe.GetType(), 
		INTEGRATION_ORDER, INTEGRATION_TYPE); 
	elmat.SetSize(trial_fe.GetNumNodes(), test_fe.GetNumNodes()); 
	elmat = 0.; 
	ElTrans& trans = trial_fe.GetTrans(); 
	double c = 1.; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		trial_fe.CalcShape(quad->X(n), _shape);
		test_fe.CalcShape(quad->X(n), _shape2); 
		_shape.OuterProduct(_shape2, _op); 
		if (_c) c = _c->Eval(trans, quad->X(n)); 
		_op *= c * quad->Weight(n) * trans.Determinant(); 
		elmat += _op; 
	}
}

void MassLumpingIntegrator::Assemble(Element& el, Matrix& elmat) {
	Quadrature* quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	ElTrans& trans = el.GetTrans(); 
	elmat.SetSize(el.GetNumNodes()); 
	elmat = 0; 
	double c = 1; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		el.CalcShape(quad->X(n), _shape);
		_shape.OuterProduct(_shape, _op);
		if (_c) c = _c->Eval(trans, quad->X(n)); 
		_op *= c * quad->Weight(n) * trans.Determinant(); 
		elmat += _op;
	}
	for (int i=0; i<elmat.Height(); i++) {
		for (int j=0; j<elmat.Width(); j++) {
			if (i!=j) {
				elmat(i,i) += elmat(i,j); 
				elmat(i,j) = 0; 
			}
		}
	}
	CHECK(elmat.IsDiagonal()); 
}

void VectorMassIntegrator::Assemble(Element& el, Matrix& elmat) {
	Quadrature* quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	elmat.SetSize(el.GetDim() * el.GetNumNodes()); 
	elmat = 0.; 
	ElTrans& trans = el.GetTrans(); 
	double c = 1.; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		el.CalcShape(quad->X(n), _shape); 
		_shape.OuterProduct(_shape, _op); 
		if (_c) c = _c->Eval(trans, quad->X(n)); 
		_op *= c * quad->Weight(n) * trans.Determinant(); 
		for (int i=0; i<el.GetDim(); i++) {
			elmat.AddMatrix(_op, i*el.GetNumNodes(), i*el.GetNumNodes()); 
		}
	}
}

void ConvectionIntegrator::Assemble(Element& el, Matrix& elmat) {
	Quadrature* quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	ElTrans& trans = el.GetTrans(); 
	elmat.SetSize(el.GetNumNodes()); 
	_v.SetSize(el.GetMeshDim()); 
	elmat = 0.; 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		el.CalcShape(quad->X(n), _shape); 
		el.CalcPhysGradShape(trans, _pgshape); 
		_vc->Eval(trans, quad->X(n), _v); 
		_v.Mult(_pgshape, _tmp); 
		_tmp *= -quad->Weight(n) * trans.Determinant(); 
		_tmp.OuterProduct(_shape, _op); 
		elmat += _op; 
	}
}

void VectorDivergenceIntegrator::MixedAssemble(Element& trial_fe, Element& test_fe, 
	Matrix& elmat) {
	Quadrature* quad = QRules.Get(trial_fe.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	elmat.SetSize(test_fe.GetNumNodes(), trial_fe.GetDim()*trial_fe.GetNumNodes()); 
	elmat = 0.; 
	ElTrans& trans = trial_fe.GetTrans(); 
	for (int n=0; n<quad->NumPoints(); n++) {
		trans.SetX(quad->X(n)); 
		test_fe.CalcShape(quad->X(n), _shape); 
		trial_fe.CalcPhysGradShape(trans, _gshape); 
		_gshape.GradToDiv(_divshape); 
		_shape *= trans.Determinant() * quad->Weight(n); 
		_shape.OuterProduct(_divshape, _op); 
		elmat += _op; 
	}
}

void UpwindFaceIntegrator::AssembleFaceMatrix(Element* e, Element* ep, 
	FaceTransformations& fts, Matrix& elmat) {
	ElTrans& etrans = e->GetTrans(); 

	Quadrature* quad = QRules.Get(fts.face->GetEl().GetType(), 
		INTEGRATION_ORDER, INTEGRATION_TYPE); 
	int size = (ep) ? e->GetNumNodes()+ep->GetNumNodes() : e->GetNumNodes(); 
	elmat.SetSize(size); 
	elmat = 0.; 
	if (ep) {
		_u.SetSize(ep->GetNumNodes()); 
	}
	_d.SetSize(e->GetNumNodes()); 
	_u = 0.; 
	_d = 0.; 
	Vector v(e->GetMeshDim()); 
	Vector nor(e->GetMeshDim()); 
	for (int i=0; i<quad->NumPoints(); i++) {
		Point ip = quad->X(i); 

		// transformations 
		fts.e_iptrans->SetX(ip); 
		Point ip_1 = fts.e_iptrans->GetPhysX(); // transform 1D ip to 2D point at boundary
		etrans.SetX(ip_1); 

		// compute normal + dot product with vector coefficient 
		fts.face->SetX(ip); 
		CalcNormal(fts.face->Jacobian(), nor); 
		_vc->Eval(etrans, ip_1, v); 
		double dot = v * nor; 
		double alpha = .5*(fabs(dot) + dot); 
		double beta = .5*(dot - fabs(dot)); 

		// downwind contribution 
		e->CalcShape(ip_1, _shape); 
		_shape.OuterProduct(_shape, _d); 
		_d *= alpha * quad->Weight(i) * fts.face->Weight(); 

		// upwind contribution 
		if (ep) {
			fts.ep_iptrans->SetX(ip); 			
			Point ip_2 = fts.ep_iptrans->GetPhysX(); 

			ep->CalcShape(ip_2, _pshape);
			_shape.OuterProduct(_pshape, _u); 
			_u *= beta * quad->Weight(i) * fts.face->Weight();  
			elmat.AddMatrix(_u, 0, e->GetNumNodes());
		}

		// add to elmat 
		elmat.AddMatrix(_d, 0, 0); 
	}
}

} // end namespace fem 