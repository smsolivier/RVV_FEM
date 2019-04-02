#include "SquareMesh.hpp"

#define FLAT(a,b) b + (Nx+1)*(a)

namespace fem 
{

SquareMesh::SquareMesh(int Nx, int Ny, Point low, Point high, Array<int> bcs) {
	_Nx = Nx; 
	_Ny = Ny; 
	_low = low; 
	_high = high; 
	int nNodes = (Nx+1) * (Ny+1); 
	int nEl = Nx * Ny; 
	_nodes.Resize(nNodes); 
	double dx = (high[0] - low[0])/(double)Nx; 
	double dy = (high[1] - low[1])/(double)Ny; 

	_dim = 2; 

	// build nodes array 
	for (int i=0; i<Ny+1; i++) {
		for (int j=0; j<Nx+1; j++) {
			MeshNode node; 
			node.id = FLAT(i,j); 
			node.x[0] = j*dx + low[0]; 
			node.x[1] = i*dy + low[1];  
			node.x[2] = 0; 
			if (j==0) node.bc = bcs[0]; 
			else if (i==Nx) node.bc = bcs[1]; 
			else if (j==Ny) node.bc = bcs[2]; 
			else if (i==0) node.bc = bcs[3]; 
			else node.bc = INTERIOR; 
			_nodes[node.id] = node;  
		}
	}

	// build elements 
	for (int i=0; i<Ny; i++) {
		for (int j=0; j<Nx; j++) {
			MeshEl el; 
			el.SetType(QUAD); 
			el.AddNode(FLAT(i,j)); 
			el.AddNode(FLAT(i,j+1)); 
			el.AddNode(FLAT(i+1,j+1)); 
			el.AddNode(FLAT(i+1,j)); 

			el.SetID(_el.GetSize()); 
			_el.Append(el); 
		}
	}

	// build boundary elements 
	for (int e=0; e<_el.GetSize(); e++) {
		MeshEl& el = _el[e]; 
		for (int n=0; n<el.GetNumNodes(); n++) {
			if (_nodes[el[n]].bc != INTERIOR) {
				int next = (n+1)%el.GetNumNodes(); 
				MeshEl nel; 
				nel.SetType(LINE); 
				nel.AddNode(el[n]); 
				nel.AddNode(el[next]); 
				nel.SetID(_bel.GetSize()); 
				_bel.Append(nel); 
			}
		}
	}

	FindNeighbors(); 
}

} // end namespace fem 