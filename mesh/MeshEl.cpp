#include "MeshEl.hpp" 

using namespace std; 
namespace fem 
{

bool in(int x, const MeshEl& el) {
	for (int i=0; i<el.GetNumNodes(); i++) {
		if (x == el[i]) return true; 
	}

	return false; 
}

int MeshEl::operator[](int i) const {
	return _nodes[i]; 
}

int MeshEl::GetTag(int i) const {
	return _tags[i]; 
}

bool MeshEl::IsNeighbor(const MeshEl& el) const {
	if (el.GetID() == (*this).GetID()) return false; 
	for (int i=0; i<_nodes.GetSize(); i++) {
		int next = (i+1)%_nodes.GetSize(); 

		if (in((*this)[i], el) && in((*this)[next], el)) {
			return true; 
		}
	}
	return false; 
}

bool MeshEl::IsFaceNeighbor(int face, const MeshEl& el) const {
	if (el.GetID() == (*this).GetID()) return false; 

	int next = (face+1)%_nodes.GetSize(); 

	if (in((*this)[face], el) && in((*this)[next], el)) {
		return true; 
	}
	return false; 
}

void MeshEl::SetNeighbors(Array<MeshNode>& nodes) {
	for (int j=0; j<GetNumNodes(); j++) {
		int k = (j+1)%GetNumNodes(); 
		MeshNode& nj = nodes[(*this)[j]]; 
		MeshNode& nk = nodes[(*this)[k]]; 
		Array<int> intersection; 
		nj.elements.Intersection(nk.elements, intersection); 
		if (intersection.GetSize() == 2) {
			if (intersection[0]==GetID()) {
				_neighbors.Append(intersection[1]); 
			} else {
				_neighbors.Append(intersection[0]); 
			}
		} else {
			_neighbors.Append(-1); 
		}
	}
}

} // end namespace fem 