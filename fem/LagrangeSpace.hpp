#pragma once 

#include "General.hpp"
#include "FESpace.hpp"
#include "Element.hpp"
#include "Polynomial.hpp"

#ifdef USE_RISCV 
extern "C" void CalcShape_RV(int Nn, int Nc, 
	const double* coef_x, const double* coef_y, double* shape, double* x); 
#endif

namespace fem 
{

/// collection of lagrange finite elements 
class LagrangeSpace : public FESpace {
public:
	/// constructor 
	LagrangeSpace(const Mesh& mesh, int order, int vdim=1); 
private:
}; 

/// represent a lagrange line element 
class LagrangeLine : public Element {
public:
	/// constructor 
	LagrangeLine(Array<MeshNode> node_list, int order=1, int mdim=2); 

	/// evaluate basis functions 
	void CalcShape(Point x, Vector& shape) const; 
	/// evaluate gradient of basis functions 
	void CalcGradShape(Point x, Matrix& gradshape) const; 
private:
	/// 1D polynomials 
	Array<Poly1D> _p; 
	/// gradients of 1D polynomials 
	Array<Poly1D> _dp; 
}; 

/// represent a lagrange quadrilateral finite element 
class LagrangeQuad : public Element {
public:
	/// constructor 
	LagrangeQuad(Array<MeshNode> node_list, int order=1, int mdim=2); 

	/// evaluate basis functions 
	void CalcShape(Point x, 
		Vector& shape) const; 
	/// evaluate gradient of basis functions 
	void CalcGradShape(Point x, 
		Matrix& gradshape) const; 
	/// evaluate gradient of basis functions in physical space 
	void CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const; 
private:
	/// 1D polynomials for x direction 
	Array<Poly1D> _p;
	/// 1D polynomials for gradient of basis function 
	Array<Poly1D> _dp; 
	/// products of 1D polynomials 
	Array<PolyProduct> _pp; 
	/// gradients 
	Array<Array<PolyProduct>> _dpp; 
	/// batched version of 1d x shape functions 
	Array<double> _shapex; 
	/// batched version of 1d y shape functions 
	Array<double> _shapey; 
	/// batched derivative of x shape functions 
	Array<double> _dshapex; 
	/// batched derivative of y shape functions 
	Array<double> _dshapey; 
}; 

/// represent a lagrange triangle finite element 
class LagrangeTri : public Element {
public:
	/// constructor 
	LagrangeTri(Array<MeshNode> node_list, int order=1, int mdim=2); 

	/// evaluate basis functions 
	void CalcShape(Point x, Vector& shape) const; 
	/// evaluate gradient of basis functions (in reference space) 
	void CalcGradShape(Point x, Matrix& gradshape) const; 
	/// evaluate gradients in physical space 
	void CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const; 
private:
	/// 1d polynomials for basis functions 
	Array<Poly1D> _p; 
	/// 1d polynomials for gradient of basis functions 
	Array<Poly1D> _dp; 
}; 

/// represent a 3D hexahedral element (cube) 
class LagrangeHex : public Element {
public: 
	/// constructor 
	LagrangeHex(Array<MeshNode> node_list, int order=1, int mdim=3); 

	/// evaluate basis functions 
	void CalcShape(Point x, Vector& shape) const; 
	/// evaluate gradient of basis functions 
	void CalcGradShape(Point x, Matrix& gradshape) const; 
	/// evaluate gradient in physical space 
	void CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const; 
private:
	/// array of 1D shape functions 
	Array<Poly1D> _p; 
	/// array of derivatives of 1D shape functions 
	Array<Poly1D> _dp; 
	/// polynomial products to create 3D shape functions 
	Array<PolyProduct> _pp; 
	/// polynomial products for 3D gradients 
	Array<Array<PolyProduct>> _dpp; 
}; 

}