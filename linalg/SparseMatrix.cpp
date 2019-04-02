#include "SparseMatrix.hpp"
#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std; 

namespace fem 
{

SparseMatrix::SparseMatrix() : Operator() {
	_zero = 0.; 
	_nnz = 0; 
}

SparseMatrix::SparseMatrix(int M, int N) : Operator(M, N) {
	_zero = 0.0; 

	_rowIndex.resize(_m); 
	_data.resize(_m); 

	_nnz = 0; 
}

void SparseMatrix::Resize(int M, int N) {
	Operator::Resize(M, N); 
	_rowIndex.resize(_m); 
	_data.resize(_m); 
}

SparseMatrix::SparseMatrix(const SparseMatrix& sp) {
	_zero = 0.; 
	_m = sp.Height(); 
	_n = sp.Width(); 
	_rowIndex.resize(_m); 
	_data.resize(_m); 
}

void SparseMatrix::Mult(const Vector& x, Vector& b) const {
	CHECK(x.GetSize() == _n); 

	if (b.GetSize() != _m) {
		b.SetSize(_m); 
		for (int i=0; i<_m; i++) {
			b[i] = 0.; 
		}		
	}

	// loop over rows 
	#pragma omp parallel for 
	for (int i=0; i<_m; i++) {
		// loop over non-zero entries in each row 
		for (int j=0; j<_rowIndex[i].size(); j++) {
			int col = _rowIndex[i][j]; 
			CHECK(col < b.GetSize()); 
			CHECK(col < x.GetSize()); 
			b[i] += _data[i][j] * x[col]; 
		}
	}
}

void SparseMatrix::operator*=(double val) {

	#pragma omp parallel for 
	for (int i=0; i<_data.size(); i++) {
		for (int j=0; j<_data[i].size(); j++) {
			_data[i][j] *= val; 
		}
	}
}

double& SparseMatrix::operator()(int row, int col) {
	CHECKMSG(row < _m, "row = " << row << ", height = " << _m); 
	CHECKMSG(col < _n, "col = " << col << ", width = " << _n); 

	// search over the columns of row for col 
	for (int i=0; i<_rowIndex[row].size(); i++) {
		// if found return 
		if (_rowIndex[row][i] == col) return _data[row][i]; 
	}

	// if not found 
	// create a new entry and return the reference to it 
	int index = _rowIndex[row].size(); // save the new index 
	_rowIndex[row].push_back(col); // append to end (order doesn't matter) 
	_data[row].push_back(_zero); // append zero 
	_nnz++; // update number of non-zeros 

	return _data[row][index]; 
}

double SparseMatrix::operator()(int row, int col) const {
	return At(row, col); 
}

double SparseMatrix::At(int row, int col) const {
	CHECK(row < _m); 
	CHECK(col < _n); 

	for (int i=0; i<_rowIndex[row].size(); i++) {
		if (_rowIndex[row][i] == col) return _data[row][i]; 
	}

	// if not found return zero 
	return _zero; 
}

void SparseMatrix::ClearNonZeros(double tol) {
	for (int i=0; i<_m; i++) {
		vector<int> new_row; 
		vector<double> new_data; 
		for (int j=0; j<_rowIndex[i].size(); j++) {
			if (abs(_data[i][j]) > tol) {
				new_row.push_back(_rowIndex[i][j]); 
				new_data.push_back(_data[i][j]); 
			}
		}

		_rowIndex[i] = new_row; 
		_data[i] = new_data; 
	}
}

void SparseMatrix::EliminateRowIntoRHS(int rc, Vector& rhs, double val) {
	CHECK(rc < Height()); 

	// subtract column rc from rhs 
	{
		for (int i=0; i<Height(); i++) {
			if (val != 0) {
				rhs[i] -= At(i,rc)*val; 			
			}
			Get(i,rc) = 0.; 
		}
	}

	rhs[rc] = val; 

	// set row rc to only one on the diagonal 
	{
		for (int i=0; i<_rowIndex[rc].size(); i++) {
			if (_rowIndex[rc][i] == rc) _data[rc][i] = 1.; 
			else _data[rc][i] = 0.; 
		}
	}
}

void SparseMatrix::Print(ostream& stream) const {
	for (int i=0; i<_m; i++) {
		for (int j=0; j<_n; j++) {
			stream << At(i, j) << " "; 
		}
		stream << endl; 
	}
}

void SparseMatrix::Sparsity(ostream& stream) const {
	stream << "SparseMatrix info:\n\tnumber of unknowns = " << _m << endl; 
	stream << "\tnumber of non-zeros = " << _nnz << endl; 
	stream << "\tpercent non-zero = " << (double)_nnz/(_m*_n) << endl; 
}

bool SparseMatrix::IsSymmetric() const {
	for (int i=0; i<_m; i++) {
		for (int j=0; j<_n; j++) {
			if (!EQUAL(At(i,j), At(j,i))) return false; 
		}
	}
	return true; 
}

void SparseMatrix::Transpose(SparseMatrix& T) const {
	T.Resize(_n, _m); 
	for (int i=0; i<_m; i++) {
		for (int j=0; j<_rowIndex[i].size(); j++) {
			T(_rowIndex[i][j], i) = _data[i][j]; 
		}
	}
}

void SparseMatrix::GetSLUFormat(double* val, int* row, int* colptr) const {
	colptr[0] = 0; 

	int index = 0; 
	for (int i=0; i<_n; i++) {
		int nvals = 0; 
		for (int j=0; j<_m; j++) {
			double el = At(i,j); 
			if (el != _zero) {
				CHECK(index < _nnz); 
				row[index] = j; 
				val[index] = el; 
				index++; 
				nvals++; 
			}
		}
		colptr[i+1] = colptr[i] + nvals; 
	}
}

#ifdef USE_EIGEN
void SparseMatrix::GetEigenFormat(Eigen::SparseMatrix<double>& eigen) const {
	eigen.resize(_m, _n); 

	for (int i=0; i<_m; i++) {
		for (int j=0; j<_rowIndex[i].size(); j++) {
			eigen.insert(i, _rowIndex[i][j]) = _data[i][j]; 
		}
	}
}
#endif

} // end namespace fem 