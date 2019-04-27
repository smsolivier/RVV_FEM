#include "FEMatrix.hpp"
#include "Opt.hpp"

using namespace std; 
namespace fem 
{

FEMatrix::FEMatrix(const FESpace* space) : Operator(space->GetVSize()) {
	_space = space; 
	int Ne = _space->GetNumElements(); 
	_data.Resize(Ne); 
	for (int i=0; i<Ne; i++) {
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
	int Ne = _space->GetNumElements(); 
#ifdef RV_MVOUTER 
	if (_mats.GetSize()==0) ERROR("must call ConvertToBatch first"); 
	int height = _data[0]->Height(); 
	Vector ball(height*Ne); 
	MVOuterC_RV(height, Ne, _mats.GetData(), 
		_vdofs.GetData(), x.GetData(), ball.GetData()); 
	BatchAdd_RV(height, Ne, _vdofs.GetData(), ball.GetData(), b.GetData()); 
#else
	Array<int> vdofs; 
	Vector elvec; 
	Vector prod; 
	for (int e=0; e<Ne; e++) {
		vdofs = _space->GetVDofs(e);  
		x.GetFromDofs(vdofs, elvec); 
		_data[e]->Mult(elvec, prod); 
		b.AddFromDofs(vdofs, prod); 
	}
#endif
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

void FEMatrix::ConvertToBatch() {
#ifdef RV_MVOUTER
	int Ne = _space->GetNumElements(); 
	int N = _data[0]->Height(); 
	_mats.Resize(N*N*Ne); 
	_vdofs.Resize(N*Ne); 
	Array<int> vdofs; 

	for (int e=0; e<Ne; e++) {
		vdofs = _space->GetVDofs(e); 
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				_mats[N*N*e + N*i + j] = _data[e]->operator()(i,j); 
			}
			_vdofs[N*e+i] = vdofs[i]; 
		}
	}
#endif
}

} // end namespace fem 