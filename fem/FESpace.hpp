#pragma once 

#include "General.hpp" 
#include "Mesh.hpp" 
#include "Node.hpp"
#include "Element.hpp"
#include "Array.hpp"

namespace fem 
{

/// represent a finite element grid 
class FESpace {
public:
	/// constructor 
	/** \param mesh mesh to build FEM grid on 
		\param order FEM order 
		\param dim treat as vector for dim > 1
	*/ 
	FESpace(const Mesh& mesh, int order, int vdim=1); 
	/// destructor 
	~FESpace(); 
	/// return the number of unknowns 
	int GetNumNodes() const {return _nodes.GetSize(); }
	/// return the total number of unknowns = _vdim * _nodes.size() 
	int GetVSize() const {return _vdim * _nodes.GetSize(); }
	/// get the global ids correspond to an element 
	/** \param[in] el element index 
		\param[out] vdofs Array of global indices (includes vector unknowns)
	*/ 
	void GetVDofs(int el, Array<int>& vdofs) const; 
	/// return the number of boundary nodes 
	int GetNBN() const {return _bnodes.GetSize(); }
	/// return the number of elements 
	int GetNumElements() const {return _el.GetSize(); }
	/// return the number of boundary elements 
	int GetNBE() const {return _bel.GetSize(); }
	/// const access to elements 
	Element& GetEl(int ind) const; 
	/// const access to boundary element ind 
	Element& GetBoundaryEl(int ind) const; 
	/// const access to nodes 
	const Node& GetNode(int ind) const; 
	/// const access to boundary nodes 
	const Node& GetBoundaryNode(int ind) const;
	/// returns all transformations needed for face integration 
	/** \param[in] e element id 
		\param[in] ep element id for neighboring element 
		\param[out] fts object with all face transformations 
	*/ 
	void GetFaceTransformations(int e, int ep, FaceTransformations& fts) const; 
	/// return the finite element order 
	int GetOrder() const {return _order; }
	/// print mesh info 
	void PrintMeshInfo(std::ostream& out=std::cout) const; 
	/// return the number of unique element tags 
	int GetNumTags() const {return _unique_tags.GetSize(); }
	/// return the ith unique tag 
	int GetTag(int i) const {
		return _unique_tags[i]; 
	}
	/// return the vector dimension 
	int GetVDim() const {return _vdim; }
protected: 
	/// store reference to mesh 
	const Mesh& _mesh; 
	/// fem order 
	int _order; 
	/// list of FEM nodes 
	Array<Node> _nodes; 
	/// list of boundary nodes 
	Array<Node> _bnodes; 
	/// list of elements 
	Array<Element*> _el; 
	/// list of elements with boundary nodes in them 
	Array<Element*> _bel; 
	/// dimensionality of problem
	int _dim; 
	/// vector dimension 
	int _vdim; 
	/// number of unique tags 
	Array<int> _unique_tags; 
}; 

} // end namespace fem 