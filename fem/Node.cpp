#include "Node.hpp" 

using namespace std; 

namespace fem 
{

Node::Node(Point x, int refID, int globalID, int bc) 
	: _x(x), _rid(refID), _gid(globalID), _bc(bc) { 
}

bool Node::CheckEqual(const Node& node) const {
	CH_TIMERS("check node equal"); 
	for (int d=0; d<DIM; d++) {
		if (!EQUAL(node.GetX()[d], _x[d])) return false; 
	}
	return true; 
}

ostream& Node::PrintX(ostream& stream) const {
	for (int i=0; i<DIM; i++) {
		stream << _x[i] << " "; 
	}
	return stream << endl; 
}

void Node::Print(ostream& out) const {
	out << "Node " << _gid << " information:" << endl; 
	out << "\tPosition: "; 
	for (int d=0; d<DIM; d++) {
		out << _x[d] << " "; 
	}
	out << "\n\tBoundary Condition = " << _bc << endl; 
	out << "\tReference ID = " << _rid << endl; 
}

} // end namespace fem 