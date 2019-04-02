#include "LHS.hpp"

using namespace std; 

namespace fem 
{

LHS::LHS(const FESpace* space, Quadrature* gq) 
	: SparseMatrix(space->GetVSize()) {

	_space = space; 
	_gq = gq; 
	_nel = _space->GetNumElements(); 
	_nnodes = _space->GetNumNodes(); 
}

void LHS::AddIntegrator(BilinearIntegrator* integ) {
	Array<int> vdofs; 
	for (int n=0; n<_nel; n++) {
		Element& el = _space->GetEl(n); 
		Matrix local(el.GetNumNodes()); 
		integ->Assemble(el, local); 
		_space->GetVDofs(n, vdofs); 
		for (int i=0; i<vdofs.GetSize(); i++) {
			for (int j=0; j<vdofs.GetSize(); j++) {
				Get(vdofs[i], vdofs[j]) += local(i,j); 
			}
		}
	}
	delete integ; 
}

void LHS::AddFaceIntegrator(BilinearIntegrator* integ) {
	Matrix elmat; 
	for (int n=0; n<_nel; n++) {
		Element& el = _space->GetEl(n); 
		for (int f=0; f<FacesPerEl(el.GetType()); f++) {
			FaceTransformations fts; 
			_space->GetFaceTransformations(n, f, fts); 
			if (el.GetNeighbor(f) >= 0) {
				Element& neighb = _space->GetEl(el.GetNeighbor(f)); 
				integ->AssembleFaceMatrix(&el, &neighb, fts, elmat); 
				Array<int> vdofs_e, vdofs_ep; 
				_space->GetVDofs(n, vdofs_e); 
				_space->GetVDofs(fts.ep_id, vdofs_ep); 
				Array<int> vdofs; 
				vdofs.Append(vdofs_e);
				vdofs.Append(vdofs_ep); 
				for (int i=0; i<vdofs.GetSize(); i++) {
					for (int j=0; j<vdofs.GetSize(); j++) {
						if (elmat(i,j) != 0.) 
							Get(vdofs[i], vdofs[j]) += elmat(i,j); 
					}
				}
			} else {
				integ->AssembleFaceMatrix(&el, NULL, fts, elmat); 
				Array<int> vdofs; 
				_space->GetVDofs(n, vdofs); 
				for (int i=0; i<vdofs.GetSize(); i++) {
					for (int j=0; j<vdofs.GetSize(); j++) {
						if (elmat(i,j) != 0.) 
							Get(vdofs[i], vdofs[j]) += elmat(i,j); 
					}
				}
			}
		}
	}
	delete integ; 
}

void LHS::ApplyDirichletBoundary(RHS& rhs, double val) {

	for (int i=0; i<_space->GetNBN(); i++) {
		const Node& node = _space->GetBoundaryNode(i); 
		if (node.GetBC() == DIRICHLET) {
			EliminateRowIntoRHS(node.GetGlobalID(), rhs, val); 
		}
	}

	// for (int i=0; i<_space->GetNumNodes(); i++) {
	// 	const Node& node = _space->GetNode(i); 
	// 	if (node.GetBC() == DIRICHLET) {
	// 		EliminateRowIntoRHS(node.GetGlobalID(), rhs, val); 
	// 	}
	// }
	ClearNonZeros(); 
}

MixedLHS::MixedLHS(const FESpace* trial, const FESpace* test, 
	Quadrature* quad) : SparseMatrix(test->GetVSize(), trial->GetVSize()) {
	_quad = quad; 
	_trial = trial; 
	_test = test; 

	_nel = _trial->GetNumElements(); 
	CHECK(_nel == _test->GetNumElements()); 
}

void MixedLHS::AddIntegrator(BilinearIntegrator* integ) {
	Array<int> vdofs_test, vdofs_trial; 
	for (int n=0; n<_nel; n++) {
		Element& test_fe = _test->GetEl(n); 
		Element& trial_fe = _trial->GetEl(n); 
		Matrix local; 
		integ->MixedAssemble(trial_fe, test_fe, local); 
		_test->GetVDofs(n, vdofs_test); 
		_trial->GetVDofs(n, vdofs_trial); 
		for (int i=0; i<vdofs_test.GetSize(); i++) {
			for (int j=0; j<vdofs_trial.GetSize(); j++) {
				Get(vdofs_test[i], vdofs_trial[j]) += local(i,j); 
			}
		}
	}
}

} // end namespace fem 