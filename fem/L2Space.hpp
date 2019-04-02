#pragma once 

#include "General.hpp"
#include "FESpace.hpp"
#include "Element.hpp"
#include "Polynomial.hpp" 
#include "Array.hpp"

namespace fem 
{

/// collection of L2 (DG) finite elements 
class L2Space : public FESpace {
public:
	/// constructor 
	L2Space(const Mesh& mesh, int order); 
}; 

/// represent an L2 segment finite element 
class L2Segment : public Element {
public:
	/// constructor 
	L2Segment(Array<MeshNode> list, int order, int mdim=1); 

	/// evaluate basis functions 
	void CalcShape(Point x, Vector& shape) const; 
	/// evaluate derivative of basis functions 
	void CalcGradShape(Point x, Matrix& gradshape) const; 
private:
	/// 1d polynomials for basis function evaluation 
	Array<Poly1D> _p; 
	/// derivative of 1d polynomials 
	Array<Poly1D> _dp; 
}; 

/// represent an L2 quadrilaterial finite element 
class L2Quad : public Element {
public:
	/// constructor 
	L2Quad(Array<MeshNode> node_list, int order, int mdim=2);

	/// evaluate basis functions 
	void CalcShape(Point x, Vector& shape) const; 
	/// evaluate gradient of basis functions 
	void CalcGradShape(Point x, Matrix& gradshape) const;
private:
	/// 1d polynomials for basis function evaluation 
	Array<Poly1D> _p; 
	/// derivative of 1D polynomials 
	Array<Poly1D> _dp; 
	Array<PolyProduct> _pp; 
	Array<Array<PolyProduct>> _dpp; 
}; 

} // end namespace fem 