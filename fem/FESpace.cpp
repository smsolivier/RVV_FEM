#include "FESpace.hpp"
#include "Quadrature.hpp"
#ifdef USE_MPI 
#include <mpi.h> 
#endif

using namespace std; 

namespace fem 
{

FESpace::~FESpace() {
	// clean up elements
	for (int i=0; i<GetNumElements(); i++) {
		delete _el[i]; 
	}
}

void FESpace::GetVDofs(int e, Array<int>& vdofs) const {
	Element& el = GetEl(e); 
	vdofs.Resize(GetVDim() * el.GetNumNodes()); 
	for (int d=0; d<GetVDim(); d++) {
		for (int i=0; i<el.GetNumNodes(); i++) {
			// vdofs[el.GetNumNodes()*d + i] = d*GetNumNodes() + el.GetNodeGlobalID(i); 
			vdofs[el.GetNumNodes()*d+i] = GetVDim()*el.GetNodeGlobalID(i) + d; 
		}
	}
}

FESpace::FESpace(const Mesh& mesh, int order, int vdim) : _mesh(mesh), _order(order) {
	_dim = _mesh.GetDim(); 
	_vdim = vdim; 
}

Element& FESpace::GetEl(int ind) const {
	return *_el[ind]; 
}

Element& FESpace::GetBoundaryEl(int ind) const {
	return *_bel[ind]; 
}

const Node& FESpace::GetNode(int ind) const {
	return _nodes[ind]; 
}

const Node& FESpace::GetBoundaryNode(int ind) const {
	return _bnodes[ind]; 
}

void FESpace::GetFaceTransformations(int e, int f, FaceTransformations& fts) const {
	Element& el = GetEl(e); 
	fts.e_id = el.GetID(); 
	fts.ep_id = el.GetNeighbor(f); 

	// setup e transformations 
	fts.e_iptrans = &el.GetFaceRefTransFromFace(f); 
	fts.face = &el.GetFaceTransFromFace(f); 

	// if ep exists (ie not a physical boundary) 
	if (fts.ep_id >= 0) {
		Element& neighb = GetEl(fts.ep_id); 
		fts.ep_iptrans = &neighb.GetFaceRefTrans(e); 
	} else {
		fts.ep_iptrans = NULL; 
	}
}

void FESpace::PrintMeshInfo(ostream& out) const {
	double max = -1.; 
	double min = numeric_limits<double>::max(); 
	double avg = 0; 
	for (int i=0; i<_el.GetSize(); i++) {
		Quadrature* quad = QRules.Get(_el[i]->GetType(), 
			INTEGRATION_ORDER, INTEGRATION_TYPE); 

		ElTrans& trans = _el[i]->GetTrans(); 

		double a = 0; 
		for (int j=0; j<quad->NumPoints(); j++) {
			trans.SetX(quad->X(j)); 
			a += fabs(quad->Weight(j)*trans.Determinant()); 
		}

		if (a > max) max = a; 
		if (a < min) min = a; 
		avg += a; 
	}
	avg /= _el.GetSize(); 

	out << "Mesh Info:" << endl; 
	out << "\tNumber of Elements = " << GetNumElements() << endl; 
	out << "\tNumber of Nodes = " << GetNumNodes() << endl; 
	out << "\tNumber of Boundary Elements = " << GetNBE() << endl; 
	out << "\tNumber of Boundary Nodes = " << GetNBN() << endl; 
	out << "\tMax Element Area = " << max << endl; 
	out << "\tMin Element Area = " << min << endl; 
	out << "\tAverage Element Area = " << avg << endl; 
	out << "\tAverage Characteristic Length = " << sqrt(avg) << endl; 
}

} // end namespace fem 
