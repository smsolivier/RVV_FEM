#include "ElTrans.hpp"

using namespace std; 

namespace fem 
{

ElTrans::ElTrans() {
	_el_set = false; 
	_x_set = false; 
	_el = NULL; 
}

ElTrans::ElTrans(Element& el) {
	SetEl(el); 
}

void ElTrans::SetEl(Element& el) {
	_el = &el; 
	_el_set = true; 
	_dim = _el->GetDim(); 
	_mdim = _el->GetMeshDim(); 
	_J.SetSize(_dim, _mdim); 
	_Jinv.SetSize(_dim, _mdim); 
}

void ElTrans::SetX(Point x) {
	_x = x; 
	_x_set = true; 
}

Point ElTrans::GetPhysX() {
	ASSERT(_x_set && _el_set); 

	Point ret; 
	Transform(_x, ret); 
	return ret; 
}

void ElTrans::Transform(const Point& x_ref, Point& x_phys) {
	ASSERT(_el_set); 
	if (_points.Height()==0) {
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
	_el->CalcGradShape(_x, _gshape); 
	// _gshape.Mult(_el->GetNodeLocationMatrix(), _J); 
	for (int i=0; i<_dim; i++) {
		for (int j=0; j<_mdim; j++) {
			_J(i,j) = 0.; 
			for (int k=0; k<_el->GetNumNodes(); k++) {
				_J(i,j) += _gshape(i,k)*(*_el)[k].GetX()[j]; 
			}
		}
	}

	return _J; 
}

const Matrix& ElTrans::InverseJacobian() {
	ASSERT(_x_set && _el_set); 
	Jacobian(); 
	_J.Inverse(_Jinv); 
	return _Jinv; 
}

double ElTrans::Determinant() {
	Jacobian(); 
	// return abs(_J.Determinant()); 
	return _J.Determinant(); 
}

double ElTrans::Weight() {
	Jacobian(); 
	return _J.Weight(); 
}

} // end namespace fem 