#include "Matrix.hpp"

#ifdef USE_RISCV 
extern "C" void MatVec_RV(int N, int M, const double* A, const double* x, double* b); 
#endif

using namespace std; 

namespace fem 
{

Matrix::Matrix() {
	_m = 0; 
	_n = 0; 
}

Matrix::Matrix(int m, int n) {
	CHECKMSG(m > 0, "m = " << m); CHECKMSG(n > 0 || n == -1, "n = " << n); 
	_m = m; 
	_n = n; 
	if (_n == -1) _n = _m; 

	Resize(_m*_n); 
}

Matrix::Matrix(const Matrix& m) {
	_m = m.Height(); 
	_n = m.Width(); 

	Resize(_m*_n); 
	for (int i=0; i<_m*_n; i++) {
		(*this)[i] = m.GetData()[i]; 
	}
}

Matrix& Matrix::operator=(const Matrix& m) {
	CHECK(m.Height() > 0 && m.Width() > 0); 
	_m = m.Height(); 
	_n = m.Width(); 
	Resize(_m*_n); 
	for (int i=0; i<_m*_n; i++) {
		(*this)[i] = m.GetData()[i]; 
	}

	return *this; 
}

double& Matrix::operator()(int i, int j) {
	CHECKMSG(i < _m && j < _n, 
		"index = (" << i << ", " << j << "), size = " << _m << " x " << _n); 
	return (*this)[i*_n+j]; 
}

double Matrix::operator()(int i, int j) const {
	CHECK(i < _m); 
	CHECK(j < _n); 
	return (*this)[i*_n+j]; 
}

void Matrix::SetSize(int m, int n) {
	CHECKMSG(m > 0, "m = " << m); CHECKMSG(n > 0 || n==-1, "n = " << n); 
	_m = m; 
	_n = n; 
	if (_n == -1) _n = _m; 
	Resize(_m*_n); 
}

void Matrix::operator*=(double val) {
	for (int i=0; i<_m*_n; i++) {
		(*this)[i] *= val; 
	}
}

void Matrix::operator+=(const Matrix& a) {
	CHECK(a.Height()==Height() && a.Width()==Width()); 

	for (int i=0; i<_m*_n; i++) {
		(*this)[i] += a.GetData()[i]; 
	}
}

void Matrix::operator-=(const Matrix& a) {
	CHECK(a.Height()==Height() && a.Width()==Width()); 

	for (int i=0; i<_m*_n; i++) {
		(*this)[i] -= a.GetData()[i]; 
	}
}

void Matrix::AddMatrix(double a, const Matrix& mat, int row, int col) {
	for (int i=0; i<mat.Height(); i++) {
		for (int j=0; j<mat.Width(); j++) {
			(*this)(i+row, j+col) += a * mat(i,j); 
		}
	}
}

void Matrix::Inverse(Matrix& inv) const {
	if (Height()==2 && Width()==2) {
		Inverse_2x2(inv); 
	} else if (Height()==3 && Width() == 3) {
		Inverse_3x3(inv); 
	} else {
		ERROR("matrix size " << Height() << "x" << Width() << " not supported"); 
	}
}

double Matrix::Determinant() const {
	if (Height()==2 && Width()==2) {
		return (*this)(0,0)*(*this)(1,1) - (*this)(1,0)*(*this)(0,1); 
	} else if (Height()==3 && Width()==3) {
		const Matrix& m = *this; 

		return m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
			m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
			m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
	} else {
		ERROR("determinant not implemented for " << Height() << "x" << Width() << 
			" matrices"); 
	}
}

double Matrix::Weight() const {
	if (Height()==Width()) {
		return Determinant(); 
	} 

	if (Height()==1 && Width()==2) {
		return sqrt((*this)[0]*(*this)[0] + (*this)[1]*(*this)[1]); 
	} else {
		ERROR("weight not implemented for size = " << Height() << 
			" x " << Width()); 
	}
}

void Matrix::Mult(const Matrix& a, Matrix& b) const {
	CHECK(Width()==a.Height()); 
	b.SetSize(Height(), a.Width()); 
	Mult(1., a, 0., b); 
}

void Matrix::Mult(double alpha, const Matrix& B, double beta, Matrix& C) const {
	CHECK(Width() == B.Height()); 
	CHECK(C.Height() == Height()); 
	CHECK(C.Width() == B.Width()); 

#ifdef USE_LAPACK 
	char trans_a = 'N'; 
	char trans_b = 'N'; 
	int M = Height(); 
	int K = Width(); 
	int N = B.Width(); 
	int LDA = M; 
	int LDB = K; 
	int LDC = M; 

	dgemm_(&trans_a, &trans_b, &M, &N, &K, &alpha, GetData(), &LDA, 
		B.GetData(), &LDB, &beta, C.GetData(), &LDC); 
#else 
	for (int i=0; i<Height(); i++) {
		for (int j=0; j<B.Width(); j++) {
			for (int k=0; k<Width(); k++) {
				C(i,j) += beta*C(i,j) + alpha*(*this)(i,k)*B(k,j); 
			}
		}
	}
#endif
}

void Matrix::AddTransMult(const Matrix& a, Matrix& b) const {
	CHECK(Height()==a.Height()); 
	CHECK(b.Width()==Width() && b.Height()==a.Width()); 

	for (int i=0; i<Width(); i++) {
		for (int j=0; j<a.Width(); j++) {
			for (int k=0; k<Height(); k++) {
				b(i,j) += (*this)(k,i)*a(k,j); 
			}
		}
	}
}

void Matrix::Mult(const Vector& x, Vector& b) const {
	CHECK(Width()==x.GetSize()); 
	b.SetSize(Height()); 

	Mult(1., x, 0., b); 
}

void Matrix::Mult(double alpha, const Vector& x, double beta, Vector& b) const {
	CHECKMSG(x.GetSize() == Width(), "size mismatch"); 
	CHECKMSG(b.GetSize() == Height(), "size mismatch"); 

#ifdef USE_RISCV 
	MatVec_RV(Height(), Width(), GetData(), x.GetData(), b.GetData()); 
#else
	for (int i=0; i<Height(); i++) {
		b[i] = beta * b[i]; 
		for (int j=0; j<Width(); j++) {
			b[i] += alpha*(*this)(i,j) * x[j]; 
		}
	}
#endif
}

bool Matrix::IsSymmetric() const {
	for (int i=0; i<Height(); i++) {
		for (int j=0; j<Width(); j++) {
			if (!EQUAL((*this)(i,j), (*this)(j,i))) return false; 
			// if (abs((*this)(i,j) - (*this)(j,i)) > 1e-12) return false; 
		}
	}
	return true; 
}

bool Matrix::operator==(const Matrix& mat) const {
	CHECK(Height() == mat.Height() && Width() == mat.Width()); 
	for (int i=0; i<Height(); i++) {
		for (int j=0; j<Width(); j++) {
			if (!EQUAL((*this)(i,j), mat(i,j))) return false; 
		}
	}
	return true; 
}

void Matrix::GradToDiv(Vector& divshape) const {
	divshape.SetSize(Height()*Width()); 
	divshape = 0.; 
	for (int d=0; d<Height(); d++) {
		for (int i=0; i<Width(); i++) {
			divshape[d*Width()+i] = (*this)(d,i); 
		}
	}
}

void Matrix::Solve(const Vector& b, Vector& x) const {
#ifdef USE_LAPACK
	x = b; 

	int N = x.GetSize(); 
	int one = 1; 
	int ipiv[x.GetSize()]; 
	int info; 
	dgesv_(&N, &one, GetData(), &N, &ipiv[0], x.GetData(), &N, &info); 
	if (info > 0) ERROR("Lapack solve issue"); 
#else 
	GaussElim(Height(), x, b); 
#endif
}

void Matrix::Transpose(Matrix& trans) const {
	trans.SetSize(Width(), Height()); 
	for (int i=0; i<Height(); i++) {
		for (int j=0; j<Width(); j++) {
			trans(j,i) = (*this)(i,j); 
		}
	}
}

double Matrix::FrobeniusNorm() const {
	double sum = 0.; 
	for (int i=0; i<GetSize(); i++) {
		sum += pow((*this)[i],2); 
	}
	return sqrt(sum); 
}

bool Matrix::IsDiagonal() const {
	bool diag = true; 
	for (int i=0; i<Height(); i++) {
		for (int j=0; j<Width(); j++) {
			if (i!=j && !EQUAL((*this)(i,j), 0)) diag = false; 
		}
	}
	return diag; 
}

void Matrix::Print(std::ostream& out) const {
	for (int i=0; i<Height(); i++) {
		for (int j=0; j<Width(); j++) {
			out << (*this)(i,j) << " "; 
		}
		out << std::endl; 
	}
}

void Matrix::Inverse_2x2(Matrix& inv) const {
	CHECK(Height()==2 && Width()==2); 

	if (inv.Height() != 2 || inv.Width() != 2) {
		inv.SetSize(2,2); 
	}

	// double det = (*this)(0,0)*(*this)(1,1) - (*this)(0,1)*(*this)(1,0);
	double det = Determinant(); 
	ERRORCHECK(det==0, "singular matrix"); 
	det = 1./det; 
	inv(0,0) = (*this)(1,1) * det; 
	inv(1,1) = (*this)(0,0) * det; 
	inv(0,1) = (*this)(0,1) * det * -1.; 
	inv(1,0) = (*this)(1,0) * det * -1.; 
}

void Matrix::Inverse_3x3(Matrix& inv) const {
	CHECK(Height()==3 && Width()==3); 
	if (inv.Height() != 3 || inv.Width() != 3) {
		inv.SetSize(3,3); 
	}

	const Matrix& m = *this; 

	double det = Determinant(); 
	ERRORCHECK(det==0, "singular matrix"); 
	det = 1./det; 

	inv(0, 0) = (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) * det;
	inv(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * det;
	inv(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * det;
	inv(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * det;
	inv(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * det;
	inv(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * det;
	inv(2, 0) = (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1)) * det;
	inv(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * det;
	inv(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * det;
}

void CalcNormal(const Matrix& J, Vector& nor, bool normalize) {
	CHECKMSG( (J.Height()==1 && J.Width()==2) || 
		(J.Height()==2 && J.Width()==3) , "J must be 1x2 or 2x3. Height = " 
		<< J.Height() << ", Width = " << J.Width()); 
	CHECKMSG( (J.Width()== nor.GetSize()) , "size of nor must match height of J. "
		<< "Height = " << J.Height() << ", size of nor = " << nor.GetSize()); 

	const double* d = J.GetData(); 
	if (J.Width() == 2) {
		nor[0] = d[1]; 
		nor[1] = -d[0]; 
	} else {
		nor[0] = d[1]*d[5] - d[2]*d[4];
		nor[1] = d[2]*d[3] - d[0]*d[5];
		nor[2] = d[0]*d[4] - d[1]*d[3];
	}

	if (normalize) {
		nor *= 1./sqrt(nor*nor); 
	}
}

int Matrix::GaussElim(int dim, Vector& x, const Vector& b) const
{
	x = b; 
	Matrix a(*this); 
	vector<int> mGaussLastRow; 
	vector<int> mGaussLastCol; 
    int ierr = 0;
    if(mGaussLastCol.size() != dim) {
       mGaussLastCol.clear();
       mGaussLastRow.clear();
       mGaussLastCol.resize(dim, dim-1);
       mGaussLastRow.resize(dim, dim-1);
       // if(mGaussPatternOn) {
       //    setGaussPattern(dim, a);
       // }
    }

    // We first do dim reduction steps.
    for(int k = 0; k < dim; k++) {
        double aSubK = a(k,k);
        double bSubK = x[k];
        if(aSubK == 0.0) {              /* check */
            //fprintf(stderr, "gauss: a[%d][%d] = 0\n", k, k);
            //cout << "gauss: a[" << k << "][" << k << "] = 0\n";
            //exit(1);
            return -1;
        }

        for(int i = k+1; i <= mGaussLastCol[k]; i++) {
            double xtemp = a(i,k);
            if(xtemp == 0.0) continue;
            a(i,k) = 0.0;                      /* eq (4.1a) */
            xtemp /= aSubK;
            for(int j = k+1; j <= mGaussLastRow[k]; j++) {
                a(i,j) -= a(k,j)*xtemp;
            }
            x[i] -= bSubK*xtemp;              /* eq (4.1c) */
        }
    }
 
    // Now we perform dim back substitutions.
    for(int i = dim-1; i >= 0; i--) {
        x[i] /= a(i,i);
        for(int k = i - 1; k >= 0; k -= 1) {
            x[k] -= a(k,i) * x[i];
        }
    }

    // Return solution in x and b vectors
    return ierr;
}

} // end namespace fem 