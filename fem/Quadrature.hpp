#pragma once 

#include "General.hpp"
#include "Vector.hpp"
#include "Array.hpp"
#include "Point.hpp"

namespace fem 
{

enum IntegrationType {
	LEGENDRE, 
	LOBATTO 
}; 

/// abstract class for numerical quadrature 
class Quadrature {
public:
	/// constructor 
	/** \param p integration order 
		\param type integration type (legendre, lobatto) 
	*/ 
	Quadrature(int p) {_p = p; }
	/// return the ith weight 
	double Weight(int i) const {return _w[i]; }
	/// return the ith integration point 
	const Point& X(int i) const {return _x[i]; }
	/// return the number of points/weights 
	int NumPoints() const {
		CHECKMSG(_x.GetSize()>0 && _w.GetSize()>0, "quadrature rule not initialized"); 
		CHECKMSG(_x.GetSize() == _w.GetSize(), "weights and points not set right"); 
		return _x.GetSize(); 
	}
	/// return the integration order 
	int GetOrder() const {return _p; }

	/// print the rule 
	void Print(std::ostream& out=std::cout) const {
		for (int i=0; i<NumPoints(); i++) {
			out << _x[i] << ", " << _w[i] << std::endl; 
		}
	}
protected:
	/// quadrature order 
	int _p; 
	/// quadrature type 
	int _type; 
	/// integration points 
	Array<Point> _x; 
	/// integration weights 
	Array<double> _w; 
}; 

/// Gauss Legendre quadrature in 1D 
class QuadLegendre : public Quadrature {
public:
	QuadLegendre(int p); 
}; 

/// Gauss Lobatto quadrature in 1D 
class QuadLobatto : public Quadrature {
public:
	QuadLobatto(int p); 
}; 

/// combine 1D rules to form 2/3D integration rule 
class QuadTensorProduct : public Quadrature {
public:
	QuadTensorProduct(Quadrature *q, int dim); 
private:
	/// number of tensor products to do 
	int _dim; 
	/// pointer to 1D rule 
	Quadrature *_q; 
}; 

/// quadrature rule for triangles 
class QuadTri : public Quadrature {
public:
	QuadTri(int p); 
}; 

/// factory class for getting quadrature rules for a given geometry, order, and type 
class QuadRules {
public:
	/// default constructor. setup arrays of possible integrators 
	QuadRules(); 
	/// destructor 
	~QuadRules(); 
	/// get a rule 
	Quadrature* Get(int geom, int order, int type=LEGENDRE); 
private:
	/// get a 1D legendre rule 
	Quadrature* GetLegendre(int order);
	/// get a 1D lobatto rule 
	Quadrature* GetLobatto(int order);  

	/// pointers for 1, 2, 3D legendre quadrature 
	Array<Array<Quadrature*>> _leg; 
	/// pointers for 1, 2, 3D lobatto quadrature 
	Array<Array<Quadrature*>> _lob; 
	/// pointers to triangle integrators 
	Array<Quadrature*> _tri; 
}; 

/// global QuadRules 
extern QuadRules QRules; 

} // end namespace fem 