#include "ElTrans.hpp"

using namespace std; 

namespace fem 
{

ElTrans::ElTrans() {
	_el_set = false; 
	_x_set = false; 
	_el = NULL; 
	_bjac = false; 
	_bJdet = false; 
	_bJinv = false; 
}

ElTrans::ElTrans(Element& el) {
	SetEl(el); 
	_bjac = false; 
	_bJdet = false;
	_bJinv = false;  
}

void ElTrans::SetEl(Element& el) {
	_el = &el; 
	_el_set = true; 
	_dim = _el->GetDim(); 
	_mdim = _el->GetMeshDim(); 
	_J.SetSize(_dim, _mdim); 
	_Jinv.SetSize(_dim, _mdim); 
	_bjac = false;
	_bJdet = false; 
	_bJinv = false;  
}

void ElTrans::SetX(Point x) {
	_x = x; 
	_x_set = true; 
	_bjac = false; 
	_bJdet = false; 
	_bJinv = false; 
}

Point ElTrans::GetPhysX() {
	ASSERT(_x_set && _el_set); 

	Point ret; 
	Transform(_x, ret); 
	return ret; 
}

void ElTrans::Transform(const Point& x_ref, Point& x_phys) {
	CH_TIMERS("transform x"); 
	ASSERT(_el_set); 
	if (_points.Height()==0) {
		CH_TIMERS("build points matrix"); 
		_points.SetSize(_mdim, _el->GetNumNodes()); 
		for (int i=0; i<_mdim; i++) {
			for (int j=0; j<_el->GetNumNodes(); j++) {
				_points(i,j) = _el->GetNode(j).GetX()[i]; 
			}
		}
	}

	_el->CalcShape(x_ref, _shape); 
	Vector v; 
	_points.Mult(_shape, v); 
	for (int i=0; i<v.GetSize(); i++) {
		x_phys[i] = v[i]; 
	}
}

const Matrix& ElTrans::Jacobian() {
	ASSERT(_x_set && _el_set); 
	if (!_bjac) {
		CH_TIMERS("jacobian"); 
		_el->CalcGradShape(_x, _gshape); 
		_gshape.Mult(_el->GetNodeLocationMatrix(), _J); 
		_bjac = true; 
	}
	return _J; 
}

const Matrix& ElTrans::InverseJacobian() {
	ASSERT(_x_set && _el_set); 
	if (!_bJinv) {
		CH_TIMERS("inverse jacobian"); 
		Jacobian(); 
		_J.Inverse(_Jinv); 
		_bJinv = true; 		
	}
	return _Jinv; 
}

double ElTrans::Determinant() {
	if (!_bJdet) {
		CH_TIMERS("jacobian determinant"); 
		Jacobian(); 
		_Jdet = _J.Determinant(); 
		_bJdet = true; 
	}
	return _Jdet; 
}

double ElTrans::Weight() {
	CH_TIMERS("jacobian weight"); 
	Jacobian(); 
	return _J.Weight(); 
}

} // end namespace fem 