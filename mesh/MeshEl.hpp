#pragma once 

#include "General.hpp" 
#include "Array.hpp"
#include "Point.hpp"

namespace fem 
{

/// represent a node in a mesh 
struct MeshNode {
	/// global id of node 
	int id; 
	/// position in DIM dimensions 
	Point x; 
	/// boundary condition type 
	int bc; 
	/// parallel global id of node 
	int pid; 
	/// list of processors this node is used by 
	Array<int> procs; 
	/// list of elements this node belongs to 
	Array<int> elements; 
}; 

/// abstract class for a mesh element 
class MeshEl {
public:
	/// default constructor 
	MeshEl() { }
	/// constructor 
	/** \param id element number in the mesh 
		\param type GMSH element type 
	*/ 
	MeshEl(int id, int type) : _id(id), _type(type) { } 
	/// access to node ids
	int operator[](int i) const;

	/// get element number for this element 
	int GetID() const {return _id; }
	/// get element type 
	int GetType() const {return _type; }
	/// get number of gmsh tags 
	int GetNumTags() const {return _tags.GetSize(); }
	/// get tag i 
	int GetTag(int i) const; 
	/// get number of nodes 
	int GetNumNodes() const {return _nodes.GetSize(); }
	/// set neighbors based on node information 
	void SetNeighbors(Array<MeshNode>& nodes); 
	/// return the neighbors array 
	const Array<int>& GetNeighbors() const {return _neighbors; }

	/// check if element is a neighbor 
	/** checks if two sequential nodes are shared */ 
	bool IsNeighbor(const MeshEl& el) const; 
	/// check if neighbor for a specific face 
	bool IsFaceNeighbor(int face, const MeshEl& el) const; 

	/// add a mesh tag 
	void AddTag(int tag) {_tags.Append(tag); }
	/// add a node 
	void AddNode(int node) {_nodes.Append(node); }
	/// set the type of the element 
	void SetType(int type) {_type = type; }
	/// set the id 
	void SetID(int id) {_id = id; }
private:
	/// element number 
	int _id; 
	/// gmsh tags 
	Array<int> _tags; 
	/// store node ids 
	Array<int> _nodes; 
	/// element type 
	int _type; 
	/// store neighbor ids 
	Array<int> _neighbors; 
}; 

} // end namespace fem 