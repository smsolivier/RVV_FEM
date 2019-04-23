#include "FEMatrix.hpp"

namespace fem 
{

FEMatrix::FEMatrix(const FESpace* space) : Operator(space->GetVSize()) {
	_space = space; 
	_data.Resize(_space->GetNumElements()); 
	for (int i=0; i<_data.GetSize(); i++) {
		Element& el = _space->GetEl(i); 
		_data[i] = new Matrix(el.GetNumNodes()); 
	}
}

FEMatrix::~FEMatrix() {
	for (int i=0; i<_data.GetSize(); i++) {
		if (_data[i]) delete _data[i]; 
	}
}

void FEMatrix::Mult(const Vector& x, Vector& b) const {
	CH_TIMERS("FEMatrix mat vec"); 
	if (b.GetSize() != Height()) b.SetSize(Height()); 

	#pragma omp parallel for 
	for (int e=0; e<_space->GetNumElements(); e++) {
		Array<int> vdofs; 
		Vector elvec; 
		const Matrix& elmat = (*this)[e]; 
		_space->GetVDofs(e, vdofs); 
		x.GetFromDofs(vdofs, elvec); 
		Vector prod; 
		elmat.Mult(elvec, prod); 
		#pragma omp critical 
		b.AddFromDofs(vdofs, prod); 
	}
}

void FEMatrix::AddIntegrator(BilinearIntegrator* integ) {
	Array<int> vdofs; 
	for (int n=0; n<_space->GetNumElements(); n++) {
		Element& el = _space->GetEl(n); 
		Matrix elmat(el.GetNumNodes()); 
		integ->Assemble(el, elmat);
		(*this)[n] += elmat;  
	}
	delete integ; 
}

void FEMatrix::ApplyDirichletBoundary(RHS& rhs, double val) {

	// eliminate into rhs 
	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& bel = _space->GetEl(e); 
		for (int n=0; n<bel.GetNumNodes(); n++) {
			if (bel[n].GetBC()==DIRICHLET) {
				Matrix& elmat = (*this)[e]; 
				for (int i=0; i<elmat.Width(); i++) {
					rhs[bel[i].GetGlobalID()] -= val * elmat(i,n); 
				}
			}
		}
	}

	// overwrite boundaries to 0 
	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& bel = _space->GetEl(e); 
		for (int n=0; n<bel.GetNumNodes(); n++) {
			if (bel[n].GetBC()==DIRICHLET) {
				Matrix& elmat = (*this)[e]; 
				for (int i=0; i<bel.GetNumNodes(); i++) {
					elmat(i,n) = 0.; 
					elmat(n,i) = 0.; 
				}
			}
		}
	}

	// put a one on the diagonal 
	Array<int> bins(_space->GetNumNodes()); 
	bins = -1;
	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& bel = _space->GetEl(e); 
		for (int n=0; n<bel.GetNumNodes(); n++) {
			if (bel[n].GetBC()==DIRICHLET) {
				Matrix& elmat = (*this)[e]; 
				if (bins[bel[n].GetGlobalID()]<0) {
					elmat(n,n) = 1.; 
					bins[bel[n].GetGlobalID()] = 1;
				}
			}
		}
	}

	// set rhs to bc value 
	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& bel = _space->GetEl(e); 
		for (int n=0; n<bel.GetNumNodes(); n++) {
			if (bel[n].GetBC()==DIRICHLET) {
				rhs[bel[n].GetGlobalID()] = val; 
			}
		}
	}
}

void FEMatrix::ConvertToSparseMatrix(SparseMatrix& spmat) const {
	spmat.Resize(_space->GetVSize()); 
	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& el = _space->GetEl(e); 
		for (int i=0; i<el.GetNumNodes(); i++) {
			for (int j=0; j<el.GetNumNodes(); j++) {
				spmat(el[i].GetGlobalID(), el[j].GetGlobalID()) += (*this)[e](i,j); 
			}
		}
	}
}

void FEMatrix::GetDiagonal(Vector& diag) const {
	diag.SetSize(Height()); 

	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& el = _space->GetEl(e); 
		const Matrix& elmat = (*this)[e]; 
		for (int i=0; i<elmat.Height(); i++) {
			diag[el[i].GetGlobalID()] += elmat(i,i); 
		}
	}
}

void FEMatrix::DiagonalPrecondition(const Vector& diag) {
	CHECK(diag.GetSize() == Height()); 
	for (int e=0; e<_space->GetNumElements(); e++) {
		Element& el = _space->GetEl(e); 
		Matrix& elmat = (*this)[e]; 
		for (int i=0; i<elmat.Height(); i++) {
			for (int j=0; j<elmat.Width(); j++) {
				elmat(i,j) /= (diag[el.GetNodeGlobalID(i)]*diag[el.GetNodeGlobalID(j)]); 
			}
		}
	}
}

void FEMatrix::operator-=(const FEMatrix& A) {
	CHECK(Height() == A.Height()); 
	for (int e=0; e<_space->GetNumElements(); e++) {
		(*this)[e] -= A[e]; 
	}
}

} // end namespace fem 