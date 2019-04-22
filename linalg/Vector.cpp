#include "Vector.hpp"

// a += b 
extern "C" void VectorAdd_RV(int N, double* a, const double* b); 
// v *= alpha 
extern "C" void VectorScale_RV(int N, double alpha, double* a); 

namespace fem 
{

Vector::Vector() {

}

Vector::Vector(int size, double val) : Array<double>(size, val) {

}

void Vector::GetFromDofs(const Array<int>& dofs, Vector& v) const {
	v.Resize(dofs.GetSize()); 
	for (int i=0; i<dofs.GetSize(); i++) {
		v[i] = (*this)[dofs[i]]; 
	}
}

void Vector::AddFromDofs(const Array<int>& dofs, Vector& v) {
	CHECKMSG(dofs.GetSize()==v.GetSize() && v.GetSize()>0, "dofs and v sizes must agree"); 
	for (int i=0; i<dofs.GetSize(); i++) {
		(*this)[dofs[i]] += v[i]; 
	}
}

Vector Vector::operator*(double val) const {
	Vector ret(GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		ret[i] = (*this)[i]*val; 
	}
	return ret; 
}

void Vector::operator*=(double val) {
	CHECKMSG(GetSize() > 0, "vector not initialized"); 
// #ifdef USE_RISCV
	// VectorScale_RV(GetSize(), val, GetData()); 
// #else
	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] *= val; 
	}
// #endif
}

void Vector::operator+=(const Vector& a) {
	CHECK(a.GetSize() == GetSize()); 

#ifdef USE_RISCV
	VectorAdd_RV(GetSize(), GetData(), a.GetData()); 
#else
	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] += a[i]; 
	}
#endif
}

void Vector::operator-=(const Vector& a) {
	CHECK(a.GetSize() == GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] -= a[i]; 
	}
}

void Vector::operator/=(const Vector& a) {
	CHECK(a.GetSize() == GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] /= a[i]; 
	}
}

void Vector::operator*=(const Vector& a) {
	CHECK(a.GetSize() == GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] *= a[i]; 
	}
}

void Vector::OuterProduct(const Vector& a, Matrix& b) const {
	if (b.Height() != GetSize() || b.Width() != a.GetSize()) {
		b.SetSize(GetSize(), a.GetSize()); 		
	}

	for (int i=0; i<GetSize(); i++) {
		for (int j=0; j<a.GetSize(); j++) {
			b(i,j) = (*this)[i] * a[j]; 
		}
	}
}

void Vector::Mult(const Matrix& a, Vector& b) const {
	CHECKMSG(GetSize()==a.Height(), 
		"vector size = " << GetSize() << 
		" must match matrix height = " << a.Height()); 

	if (b.GetSize() != a.Width()) {
		b.Resize(a.Width()); 
	}

	b = 0.; 
	for (int i=0; i<a.Width(); i++) {
		for (int j=0; j<a.Height(); j++) {
			b[i] += (*this)[j] * a(j,i); 
		}
	}
}

void Vector::Transpose() {
	Vector tmp = (*this); 
	for (int i=0; i<tmp.GetSize(); i++) {
		(*this)[i] = tmp[tmp.GetSize()-1-i]; 
	}
}

void Vector::SquareRoot() {
	for (int i=0; i<GetSize(); i++) {
		CHECKMSG((*this)[i] >= 0, "negative value in square root"); 
		(*this)[i] = sqrt((*this)[i]); 
	}
}

double Vector::Dot(const Vector& x) const {
	CHECK(x.GetSize() == GetSize()); 

	double sum = 0; 
	// #pragma omp parallel for reduction(+:sum) 
	for (int i=0; i<GetSize(); i++) {
		sum += (*this)[i]*x[i]; 
	}

	return sum; 
}

double Vector::L2Norm() const {
	double sum = 0; 

	// #pragma omp parallel for reduction(+:sum) 
	for (int i=0; i<GetSize(); i++) {
		sum += (*this)[i]*(*this)[i]; 
	}
	return sqrt(sum); 
}

double Vector::LinfNorm() const {
	double max = 0; 
	for (int i=0; i<GetSize(); i++) {
		if ((*this)[i] > max) max = (*this)[i]; 
	}
	return max; 
}

void Vector::Linspace(int N, double a, double b) {
	double dx = (b - a)/(double)(N-1); 
	Resize(N); 
	for (int i=0; i<N; i++) {
		(*this)[i] = dx*i + a; 
	}
}

void Vector::Print(std::ostream& out, std::string end) const {
	for (int i=0; i<GetSize(); i++) {
		out << (*this)[i] << end; 
	}
	out << std::endl; 
}

void Add(const Vector& a, const Vector& b, Vector& c) {
	CHECK(a.GetSize() == b.GetSize()); 
	if (c.GetSize() != a.GetSize()) c.Resize(a.GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<a.GetSize(); i++) {
		c[i] = a[i] + b[i]; 
	}
}

void Subtract(const Vector& a, const Vector& b, Vector& c) {
	CHECK(a.GetSize() == b.GetSize()); 
	if (c.GetSize() != a.GetSize()) c.Resize(a.GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<a.GetSize(); i++) {
		c[i] = a[i] - b[i]; 
	}
}

void Add(double alpha, const Vector& a, double beta, const Vector& b, Vector& c) {
	CHECK(a.GetSize() == b.GetSize()); 
	c.Resize(a.GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<a.GetSize(); i++) {
		c[i] = alpha*a[i] + beta*b[i]; 
	}
}

} // end namespace fem 