#pragma once 
#include "General.hpp"
#include "MeshEl.hpp"

namespace fem 
{

/// represent a node of a finite element 
class Node {
public:
	/// default constructor 
	Node() { }
	/// constructor 
	/** \param x position in DIM dimensions \n
		\param refID node number in the reference element \n
		\param globaID global id number on the FEM mesh 
		\param bc boundary type as defined in Mesh 
	*/ 
	Node(Point x, int refID, int globalID, int bc); 
	/// set the global id 
	void SetGlobalID(int id) {_gid = id; }
	/// check if two nodes are equal (same x) 
	bool CheckEqual(const Node& node) const; 

	/// return location of node 
	const Point& GetX() const {return _x; }
	/// return reference id 
	int GetRefID() const {return _rid; }
	/// return global id 
	int GetGlobalID() const {return _gid; }
	/// return boundary condition 
	int GetBC() const {return _bc; }
	/// returns the processors this belongs to 
	const Array<int>& GetProcs() const {return _procs; }
	/// sets the processors that this node belongs to 
	void SetProcs(Array<int>& procs) {_procs.Append(procs); }

	/// write the position to the stream 
	std::ostream& PrintX(std::ostream& stream=std::cout) const; 
	/// write node information to the stream 
	void Print(std::ostream& out=std::cout) const; 
private:
	/// location of node in physical space 
	Point _x; 
	/// node number in the reference element 
	int _rid; 
	/// node number in global mesh 
	int _gid; 
	/// boundary type 
	int _bc; 
	/// store which processors use this node 
	Array<int> _procs; 
}; 

} // end namespace fem 