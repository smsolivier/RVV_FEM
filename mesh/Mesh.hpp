#pragma once 
#include "General.hpp"
#include "MeshEl.hpp"
#include "Array.hpp"
// #include "Node.hpp"

#define NGmshTypes 7

namespace fem 
{

/// define the possible GMSH mesh element types 
enum GmshElType {
	LINE, 
	TRI, 
	QUAD, 
	TET, 
	HEX, 
	PRISM, 
	PYR
}; 

/// define boundary conditions 
enum BoundaryCondition {
	INTERIOR, 
	DIRICHLET, 
	NEUMANN
};  

/// represent a FEM mesh 
class Mesh { 
public: 
	/// default constructor 
	Mesh() { } 
	/// load the mesh file 
	/** \param name name of mesh file (relative to the mesh directory) */ 
	Mesh(std::string name, int nref=0); 
	/// refine the mesh 
	void GlobalRefine(); 

	/// const access to MeshEls 
	const MeshEl& GetElement(int index) const; 
	/// const access to boundary elements 
	const MeshEl& GetBoundaryElement(int index) const; 
	/// const access to MeshNodes 
	const MeshNode& GetNode(int index) const; 
	/// access to MeshNodes 
	MeshNode& GetNode(int index) {return _nodes[index]; }
	/// get all MeshNodes for MeshEl ind 
	void GetNodes(int ind, Array<MeshNode>& nodes) const; 

	/// write the mesh to gmsh format 
	void Write() const {ERROR("not implemented"); }

	/// return the number of nodes 
	int GetNumNodes() const {return _nodes.GetSize(); }
	/// return the number of elements 
	int GetNumElements() const {return _el.GetSize(); }
	/// return the number of boundary elements 
	int GetNumBoundaryElements() const {return _bel.GetSize(); }
	/// return the dimensionality of the mesh 
	int GetDim() const {return _dim; }
protected:
	/// find neighbors of all elements 
	void FindNeighbors(); 
	/// number of nodes per element for each Gmsh element type 
	std::array<int,NGmshTypes> _NodesPerEl; 
	/// store all nodes 
	Array<MeshNode> _nodes; 
	/// store all elements 
	Array<MeshEl> _el; 
	/// store the boundary elements (ie for quad these will be segments) 
	Array<MeshEl> _bel; 
	/// mesh version 
	double _mesh_version; 
	/// ascii type 
	int _ascii; 
	/// size of data in mesh 
	int _data_size; 
	/// dimensionality of mesh 
	int _dim; 
}; 

/// returns the dimensionality of a given shape 
int GetGeoDim(int geo); 
/// returns the number of faces for a given element shape 
int FacesPerEl(int geo); 

} // end namespace FEM 