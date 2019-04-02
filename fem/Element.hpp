#pragma once 

#include "General.hpp"
#include "Mesh.hpp"
#include "Node.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "ElTrans.hpp"
#include "Quadrature.hpp"

namespace fem 
{

class ElTrans; 

/// base class for finite elements 
class Element {
public:
	/// default constructor
	Element() {_trans = NULL; } 
	/// construct an element from the Mesh 
	/** \param node_list list of MeshNodes to build FEM nodes on 
		\param order polynomial order for this element 
		\param geo GMSH element type 
	*/ 
	Element(Array<MeshNode>& node_list, int order, int geo, int mdim); 
	/// deconstructor 
	~Element(); 
	/// set tags for material properties etc 
	void SetTag(int tag) {_tag = tag; }
	/// return the tag for this element 
	int GetTag() const {return _tag; }
	/// set the id number of this element 
	void SetID(int id) {_id = id; }
	/// return the id number of this element 
	int GetID() const {return _id; }

	/// access to FEM nodes in element 
	Node& GetNode(int index); 
	/// access to FEM nodes in element 
	Node& operator[](int index) {return GetNode(index); }
	/// const access to FEM nodes in element 
	const Node& GetNode(int index) const; 
	/// return the global id of node i 
	int GetNodeGlobalID(int i) const; 
	/// const access to FEM nodes in element 
	const Node& operator[](int index) const {return GetNode(index); }
	/// return number of FEM nodes in this element 
	int GetNumNodes() const {return _nodes.GetSize(); }
	/// return dimensionality of the element shape 
	int GetDim() const {return _dim; }
	/// return dimensionality of mesh 
	int GetMeshDim() const {return _mdim; }
	/// return the element geometry type 
	int GetType() const {return _geo; }
	/// return the element transformation 
	ElTrans& GetTrans() {return *_trans; }
	/// set neighbor for face i 
	void SetNeighbor(int i, int ID) {
		_neighbors[i] = ID; 
	}
	/// get neighbor on face i 
	int GetNeighbor(int i) const {
		return _neighbors[i]; 
	}
	/// get the node location matrix (Nnodes x dim) 
	const Matrix& GetNodeLocationMatrix(); 
	/// get face transformation that corresponds to neighbor element e
	/** \param e element number of neighboring element 
		finds which face is corresponds to element e and returns a transformation 
		that can be used to transform _dim-1 dimensional integration rules to 
		_dim dimensional integration points for use in CalcShape 
		Face transformations should be set in constructor of derived classes 
	*/ 
	ElTrans& GetFaceRefTrans(int e); 
	/// get the FaceRefTrans from the face number 
	ElTrans& GetFaceRefTransFromFace(int f); 
	/// get physical space face transformation for face shared with element e 
	ElTrans& GetFaceTrans(int e); 
	/// get FaceTrans from face number 
	ElTrans& GetFaceTransFromFace(int f); 
	/// integrate this element 
	double Integrate(const Vector& x, Quadrature* quad=NULL) const; 
	/// return the centroid of the element 
	Point Centroid() const; 
	/// return the area/volume of this element 
	double Volume(Quadrature* quad=NULL); 
	/// return the energy norm \f$ E = \int |u|^2 dV + \int \nabla u \cdot \nabla u dV \f$ 
	double EnergyNorm(const Vector& x, Quadrature* quad=NULL) const; 
	/// return \f$ \sum B(\xi_i) u_i \f$ 
	double Interpolate(const Vector& x, const Point& p) const; 
	/// interpolate the gradient 
	/** \param[in] x vector of all node values 
		\param[in] p point to interpolate to 
		\param[out] grad interpolated gradient 
	*/ 
	void InterpolateGradient(const Vector& x, const Point& p, Vector& grad) const; 

	// virtual functions 
	/// evaluate basis function at x 
	virtual void CalcShape(Point x, 
		Vector& shape) const {
		ERROR("CalcShape not defined"); 
	} 
	/// evaluate basis functions in physical space 
	virtual void CalcPhysShape(ElTrans& trans, Vector& shape) const {
		ERROR("CalcPhysShape not defined"); 
	}
	/// evaluate gradient of basis functions at x 
	virtual void CalcGradShape(Point x, 
		Matrix& gradshape) const {
		ERROR("CalcGradShape not defined"); 
	}
	/// evaluate gradient of basis functions in physical space using ElTrans 
	virtual void CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const; 

	/// print Element information 
	void Print(std::ostream& out=std::cout) const; 
protected: 
	/// find the index into _neighbors that matches element e 
	int FindNeighbor(int e) const; 

	/// id number for this element 
	int _id; 
	/// list of nodes. These nodes do not need to correlate to mesh nodes 
	Array<Node> _nodes; 
	/// list of MeshNodes that define the geometry of the Element 
	Array<MeshNode> _geo_nodes; 
	/// polynomial order of element 
	int _order; 
	/// dimensionality of element (ie the dimension of this element shape) 
	int _dim; 
	/// dimensionality of mesh (ie the total problem). Allows representations of faces in higher dimensions 
	int _mdim; 
	/// store global id of neighboring elements 
	Array<int> _neighbors; 
	/// store the geometry type 
	int _geo; 
	/// store the material property tag 
	int _tag; 
	/// store the volume/area of this element 
	double _volume; 
	/// store the ElTrans transformation object 
	ElTrans* _trans; 
	/// store reference space face transformations for each face 
	Array<Element*> _faceRefTrans; 
	/// physical space face transformations 
	Array<Element*> _faceTrans; 
	/// store the node locations in a matrix (Nnodes x dim) 
	Matrix _nodeLocMatrix; 
}; 

} // end namespace fem 