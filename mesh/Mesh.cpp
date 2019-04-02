#include "Mesh.hpp"
#include <fstream> 
#include <limits> 
// #include "Quadrature.hpp"
// #include "ElTrans.hpp"

#ifdef USE_METIS 
#include "metis.h"
#endif

#ifdef USE_MPI 
#include <mpi.h>
#endif

#define GEO_TOL 1e-5

using namespace std; 

bool in(string str1, string str2) {
	return str1.find(str2) != string::npos; 
}

namespace fem
{

bool matchingNode(const MeshNode& n1, const MeshNode& n2) {
	for (int i=0; i<DIM; i++) {
		if (abs(n1.x[i] - n2.x[i]) > GEO_TOL) return false; 
	}
	return true; 
}

Mesh::Mesh(string name, int nref) {
	ifstream meshstream(MESH_DIR + name); 
	if (!meshstream.good()) ERROR("meshfile " << name << " not read"); 

	_NodesPerEl = {2, 3, 4, 4, 8, 6, 5}; 

	string line; 
	while (getline(meshstream, line)) {
		if (in(line, "$MeshFormat")) {
			meshstream >> _mesh_version >> _ascii >> _data_size; 
		} 

		else if (in(line, "$Nodes")) {
			int nNodes; 
			meshstream >> nNodes; 
			_nodes.Resize(nNodes); 
			for (int i=0; i<nNodes; i++) {
				MeshNode node; 
				meshstream >> node.id >> node.x[0] >> node.x[1] >> node.x[2];
				node.id--; // decrement for 0 based indexing 
				_nodes[i] = node; 
			}
		}

		else if (in(line, "$Elements")) {
			int nEl; 
			meshstream >> nEl; 
			for (int i=0; i<nEl; i++) {
				int id, type, ntags; 
				meshstream >> id >> type >> ntags; 
				type--; // decrement for 0 based indexing 

				MeshEl el; 
				el.SetType(type); 

				// set dim 
				_dim = 2; 
				if (type == TET || type == HEX || type == PRISM || type == PYR) {
					_dim = 3; 
				}

				// read in tags 
				for (int j=0; j<ntags; j++) {
					int tag; 
					meshstream >> tag;
					el.AddTag(tag);  
				}

				// read in nodes 
				int nnodes = _NodesPerEl[el.GetType()]; 
				for (int j=0; j<nnodes; j++) {
					int node; 
					meshstream >> node; 
					node--;
					el.AddNode(node); 
				}

				// set boundary nodes 
				if (el.GetType() == LINE) {
					for (int j=0; j<nnodes; j++) {
						if (el.GetTag(0)==1) {
							_nodes[el[j]].bc = DIRICHLET; 
						}
						else if (el.GetTag(0)==2) {
							_nodes[el[j]].bc = NEUMANN; 
						}
						else ERROR("boundary value " 
							<< el.GetTag(0) << " not defined"); 
					}
					// build boundary element 
					el.SetID(_bel.GetSize()); 
					_bel.Append(el); 
				} 

				// add to _el and set the element id for each node 
				else {
					el.SetID(_el.GetSize()); 
					for (int j=0; j<nnodes; j++) {
						_nodes[el[j]].elements.Append(el.GetID()); 
					}
					_el.Append(el); 
				}
			}
		}
	}

	meshstream.close();

	for (int i=0; i<nref; i++) {
		GlobalRefine(); 
	}

	FindNeighbors(); 
}

void Mesh::GlobalRefine() {

	Array<MeshEl> old_els = _el; 
	vector<vector<int>> node_map(_el.GetSize()); 
	int start_nodes = _nodes.GetSize(); 
	for (int n=0; n<old_els.GetSize(); n++) {
		// access element n
		const MeshEl& old = _el[n]; 

		if (old.GetType() != QUAD) {
			WARNING("only quads supported"); 
			return; 
		}

		// construct new nodes 
		Array<MeshNode> new_nodes(5); 
		for (int i=0; i<4; i++) {
			int next = (i+1)%old.GetNumNodes(); 
			Point x; 
			for (int j=0; j<DIM; j++) {
				x[j] = .5*(_nodes[old[i]].x[j] + _nodes[old[next]].x[j]); 
			}
			new_nodes[i].x = x; 
			if (_nodes[old[i]].bc == _nodes[old[next]].bc) {
				new_nodes[i].bc = _nodes[old[i]].bc; 
			} else {
				new_nodes[i].bc = INTERIOR; 
			}
		}

		// average node positions 
		Point x; 
		for (int i=0; i<DIM; i++) {
			x[i] = .5*(_nodes[old[0]].x[i] + _nodes[old[2]].x[i]); 
		}
		new_nodes[4].x = x; 
		new_nodes[4].bc = INTERIOR; 

		// set to node id's to -1 
		for (int i=0; i<new_nodes.GetSize(); i++) {
			new_nodes[i].id = -1; 
		}

		// loop over faces 
		for (int i=0; i<old.GetNumNodes(); i++) {
			const MeshNode& node = _nodes[old[i]]; 
			// check if neighbor on face i 
			for (int j=0; j<node.elements.GetSize(); j++) {

				// if a neighbor and neighbor has already been refined 
				int el_num = node.elements[j]; 
				if (old.IsFaceNeighbor(i, old_els[el_num])) {

					// find match for new_node[i] 
					for (int k=0; k<node_map[el_num].size(); k++) {
						// if match set node id to already existing node's id 
						int node_num = node_map[el_num][k]; 
						if (matchingNode(new_nodes[i], _nodes[node_num])) {
							new_nodes[i].id = node_num; 
						}
					}
				}				
			}
		}

		// add new nodes to _nodes 
		for (int i=0; i<new_nodes.GetSize(); i++) {

			// if a redundant node leave node id as is and 
			// don't append to end of _nodes 
			if (new_nodes[i].id > 0) {

			} else {
				new_nodes[i].id = _nodes.GetSize(); 
				_nodes.Append(new_nodes[i]); 				
			}
		}

		// construct new elements 
		vector<MeshEl> new_els(4); 
		// add first two for all elements 
		for (int i=0; i<new_els.size(); i++) {
			new_els[i].AddNode(old[i]); 
			new_els[i].AddNode(new_nodes[i].id); 
		}

		// el0 remaining 
		new_els[0].AddNode(new_nodes[4].id); 
		new_els[0].AddNode(new_nodes[3].id); 

		// el1 remaining 
		new_els[1].AddNode(new_nodes[4].id); 
		new_els[1].AddNode(new_nodes[0].id); 

		// el2 remaining 
		new_els[2].AddNode(new_nodes[4].id); 
		new_els[2].AddNode(new_nodes[1].id); 

		// el3 remaining 
		new_els[3].AddNode(new_nodes[4].id); 
		new_els[3].AddNode(new_nodes[2].id);

		// copy meta info from old element 
		for (int i=0; i<new_els.size(); i++) {
			for (int j=0; j<old.GetNumTags(); j++) {
				new_els[i].AddTag(old.GetTag(j)); 
			}
			new_els[i].SetType(old.GetType()); 
		}

		// put el0 into old's location 
		new_els[0].SetID(n); 
		_el[n] = new_els[0]; 

		// append the rest to the end 
		for (int i=1; i<new_els.size(); i++) {
			new_els[i].SetID(_el.GetSize()); 
			_el.Append(new_els[i]); 
		}

		// update the node map 
		for (int i=0; i<new_nodes.GetSize(); i++) {
			node_map[n].push_back(new_nodes[i].id); 
		}
	}

	// update elements nodes belong to 
	for (int n=0; n<_el.GetSize(); n++) {
		for (int j=0; j<_el[n].GetNumNodes(); j++) {
			_nodes[_el[n][j]].elements.Clear(); 
		}
	}

	for (int n=0; n<_el.GetSize(); n++) {
		for (int j=0; j<_el[n].GetNumNodes(); j++) {
			_nodes[_el[n][j]].elements.Append(_el[n].GetID()); 
		}
	}
}

const MeshEl& Mesh::GetElement(int index) const {
	return _el[index]; 
}

const MeshEl& Mesh::GetBoundaryElement(int index) const {
	return _bel[index]; 
}

const MeshNode& Mesh::GetNode(int index) const {
	return _nodes[index]; 
}

void Mesh::GetNodes(int ind, Array<MeshNode>& nodes) const {
	const MeshEl& el = GetElement(ind); 
	nodes.Resize(el.GetNumNodes()); 
	for (int i=0; i<el.GetNumNodes(); i++) {
		nodes[i] = GetNode(el[i]); 
	}
}

void Mesh::FindNeighbors() {
	for (int i=0; i<GetNumElements(); i++) {
		_el[i].SetNeighbors(_nodes); 
	}
}

int GetGeoDim(int geo) {
	if (geo == LINE) return 1; 
	else if (geo == QUAD || geo == TRI) return 2; 
	else if (geo == TET || geo == HEX || geo == PRISM || geo == PYR) return 3; 
	else ERROR("geometry " << geo << " not defined"); 
}

int FacesPerEl(int geo) {
	if (geo == QUAD) return 4; 
	else if (geo == TRI) return 3; 
	else if (geo == LINE) return 2; 
	else if (geo == HEX) return 6; 
	else ERROR("number of faces for geo = " << geo << " not defined"); 
	return -1; 
}

} // end namespace FEM 