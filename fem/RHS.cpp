#include "RHS.hpp"
#ifdef USE_MPI 
#include <mpi.h>
#endif

using namespace std; 

namespace fem 
{

RHS::RHS(const FESpace* space, Quadrature* gq) {
	_space = space; 
	_nnodes = _space->GetNumNodes(); 
	SetSize(_space->GetVSize()); 
}

void RHS::AddIntegrator(LinearIntegrator* integ) {
	for (int n=0; n<_space->GetNumElements(); n++) {
		Element& el = _space->GetEl(n); 
		Vector local(el.GetNumNodes()); 
		integ->Assemble(el, local); 
		for (int i=0; i<el.GetNumNodes(); i++) {
			const Node& ni = el.GetNode(i); 
			this->operator[](ni.GetGlobalID()) += local[i]; 
		}
	}
	delete integ; 
}

void RHS::AddBoundaryIntegrator(LinearIntegrator* integ) {
	Vector local;
	for (int n=0; n<_space->GetNumElements(); n++) {
		Element& el = _space->GetEl(n); 
		for (int f=0; f<FacesPerEl(el.GetType()); f++) {
			if (el.GetNeighbor(f) < 0) {
				ElTrans& face = el.GetFaceTransFromFace(f);
				Element& fel = face.GetEl(); 
				integ->Assemble(fel, local); 
				(*this)[el.GetNodeGlobalID(f)] += local[0]; 
				(*this)[el.GetNodeGlobalID((f+1)%2)] += local[1]; 
			}
		}
	}
	delete integ; 
}

} // end namespace fem 