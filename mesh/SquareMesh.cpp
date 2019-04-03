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

CubeMesh::CubeMesh(Array<int> N, Point low, Point high, Array<int> bc) {
	CHECK(N.GetSize() == DIM); 
	_n = N; 
	_low = low; 
	_high = high; 
	_dim = 3; 

	int Nn = 1; 
	int Ne = 1; 
	for (int d=0; d<N.GetSize(); d++) {
		Nn *= (N[d] + 1); 
		Ne *= N[d]; 
	}

	Point delta; 
	for (int d=0; d<DIM; d++) {
		delta[d] = (high[d] - low[d])/N[d]; 
	}

	// build nodes array 
	auto flat = [N](int i, int j, int k) {return k + j*(N[0]+1) + i*(N[0]+1)*(N[1]+1); }; 
	_nodes.Resize(Nn); 
	for (int i=0; i<N[2]+1; i++) {
		for (int j=0; j<N[1]+1; j++) {
			for (int k=0; k<N[0]+1; k++) {
				MeshNode node; 
				node.id = flat(i,j,k); 
				node.x[0] = k*delta[0] + low[0]; 
				node.x[1] = j*delta[1] + low[1]; 
				node.x[2] = i*delta[2] + low[2]; 

				if (j==0) node.bc = bc[0]; 
				else if (k==N[0]) node.bc = bc[1];
				else if (j==N[1]) node.bc = bc[2]; 
				else if (k==0) node.bc = bc[3]; 
				else if (i==N[2]) node.bc = bc[4]; 
				else if (i==0) node.bc = bc[5];  
				else node.bc = INTERIOR; 
				_nodes[node.id] = node; 
			}
		}
	}

	// build elements 
	for (int i=0; i<N[2]; i++) {
		for (int j=0; j<N[1]; j++) {
			for (int k=0; k<N[0]; k++) {
				MeshEl el; 
				el.SetType(HEX); 
				el.AddNode(flat(i,j,k)); 
				el.AddNode(flat(i,j,k+1)); 
				el.AddNode(flat(i,j+1,k+1)); 
				el.AddNode(flat(i,j+1,k)); 

				el.AddNode(flat(i+1,j,k)); 
				el.AddNode(flat(i+1,j,k+1)); 
				el.AddNode(flat(i+1,j+1,k+1)); 
				el.AddNode(flat(i+1,j+1,k)); 

				el.SetID(_el.GetSize()); 
				_el.Append(el); 
			}
		}
	}
	FindNeighbors(); 
}

} // end namespace fem 