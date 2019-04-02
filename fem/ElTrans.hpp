#pragma once 

#include "General.hpp"
#include "Element.hpp"
#include "Matrix.hpp"

namespace fem 
{

class Element; 
class ElTrans; 

/// store all the information and transformations needed for a face integral 
struct FaceTransformations {
	/// element id
	int e_id; 
	/// id of neighboring element 
	int ep_id; 
	/// ip transformation for main element
	ElTrans* e_iptrans;  
	/// ip transformation for neighbor element
	ElTrans* ep_iptrans;  
	/// transformation for the face element 
	ElTrans* face; 
}; 

/// transform from reference to physical space 
class ElTrans {
public:
	/// constructor 
	ElTrans();
	/// constructor 
	ElTrans(Element& el);  
	/// set element 
	void SetEl(Element& el); 
	/// return the element 
	Element& GetEl() const {return *_el; }
	/// set evaluation point 
	void SetX(Point x); 
	/// return evaluation point 
	const Point& GetX() const {return _x; }
	/// return evaluation point in physical space 
	/** evaluate \f$ x(r_e) = \sum_i B_i(r_e) x_{p,i} \f$ */ 
	Point GetPhysX(); 
	/// transform to physical space without saving 
	void Transform(const Point& x_ref, Point& x_phys); 
	/// evaluate jacobian 
	const Matrix& Jacobian(); 
	/// evaluate inverse jacobian 
	const Matrix& InverseJacobian(); 
	/// determinant of transformation 
	double Determinant(); 
	/// weight of transformation 
	double Weight(); 
private:
	/// store element (for access to CalcGradShape) 
	Element* _el; 
	/// true if element has been set 
	bool _el_set; 
	/// evaluation point 
	Point _x; 
	/// true if evaluation point has been set 
	bool _x_set; 
	/// store grad shape 
	Matrix _gshape; 
	/// store shape 
	Vector _shape; 
	/// store jacobian for reuse 
	Matrix _J; 
	/// store inverse jacobian for reuse 
	Matrix _Jinv; 
	/// dimensionality of element shape 
	int _dim; 
	/// dimensionality of mesh 
	int _mdim; 
	/// (dim x NumNodes) matrix of point values in physical space 
	Matrix _points; 
}; 

} // end namespace fem 