#include "Vector.hpp"
#include "Opt.hpp"

namespace fem 
{

Vector::Vector() {

}

Vector::Vector(int size, double val) : Array<double>(size, val) {

}

void Vector::GetFromDofs(const Array<int>& dofs, Vector& v) const {
	CH_TIMERS("get from dofs"); 
	v.Resize(dofs.GetSize()); 
#ifdef RV_GETDOFS
	GetFromDofs_RV(dofs.GetSize(), dofs.GetData(), v.GetData(), GetData()); 
#else
	for (int i=0; i<dofs.GetSize(); i++) {
		v[i] = (*this)[dofs[i]]; 
	}
#endif
}

void Vector::AddFromDofs(const Array<int>& dofs, Vector& v) {
	CH_TIMERS("add from dofs"); 
	CHECKMSG(dofs.GetSize()==v.GetSize() && v.GetSize()>0, "dofs and v sizes must agree"); 
#ifdef RV_ADDDOFS
	AddFromDofs_RV(dofs.GetSize(), dofs.GetData(), v.GetData(), GetData()); 
#else
	for (int i=0; i<dofs.GetSize(); i++) {
		(*this)[dofs[i]] += v[i]; 
	}
#endif
}

void Vector::operator*=(double val) {
	CH_TIMERS("vector operation *="); 
	CHECKMSG(GetSize() > 0, "vector not initialized"); 
#ifdef RV_VECSCALE
	VectorScale_RV(GetSize(), GetData(), &val); 
#else
	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] *= val; 
	}
#endif
}

void Vector::operator+=(const Vector& a) {
	CH_TIMERS("vector +="); 
	CHECK(a.GetSize() == GetSize()); 
#ifdef RV_VECADD
	VectorAdd_RV(GetSize(), GetData(), a.GetData()); 
#else
	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] += a[i]; 
	}
#endif
}

void Vector::operator-=(const Vector& a) {
	CH_TIMERS("vector -="); 
	CHECK(a.GetSize() == GetSize()); 
#ifdef RV_VECSUB
	VectorSub_RV(GetSize(), GetData(), a.GetData()); 
#else

	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] -= a[i]; 
	}
#endif
}

void Vector::operator/=(const Vector& a) {
	CH_TIMERS("vector /="); 
	CHECK(a.GetSize() == GetSize()); 
#ifdef RV_VECDIV
	VectorDiv_RV(GetSize(), GetData(), a.GetData()); 
#else
	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] /= a[i]; 
	}
#endif
}

void Vector::operator*=(const Vector& a) {
	CH_TIMERS("vector vector *="); 
	CHECK(a.GetSize() == GetSize()); 
#ifdef RV_VECMUL
	VectorMul_RV(GetSize(), GetData(), a.GetData()); 
#else
	#pragma omp parallel for 
	for (int i=0; i<GetSize(); i++) {
		(*this)[i] *= a[i]; 
	}
#endif
}

void Vector::OuterProduct(const Vector& a, Matrix& b) const {
	CH_TIMERS("outer product"); 
	if (b.Height() != GetSize() || b.Width() != a.GetSize()) {
		b.SetSize(GetSize(), a.GetSize()); 		
	}
#ifdef RV_VECOP 
	VectorOP_RV(GetSize(), a.GetSize(), GetData(), a.GetData(), b.GetData()); 
#else
	for (int i=0; i<GetSize(); i++) {
		for (int j=0; j<a.GetSize(); j++) {
			b(i,j) = (*this)[i] * a[j]; 
		}
	}
#endif
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
	CH_TIMERS("vector dot product"); 
	CHECK(x.GetSize() == GetSize()); 
#ifdef RV_VECDOT
	double ret = 0; 
	VectorDot_RV(GetSize(), GetData(), x.GetData(), &ret); 
	return ret; 
#else
	double sum = 0; 
	// #pragma omp parallel for reduction(+:sum) 
	for (int i=0; i<GetSize(); i++) {
		sum += this->operator[](i)*x[i]; 
	}
	return sum; 
#endif
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
	CH_TIMERS("vvadd"); 
	CHECK(a.GetSize() == b.GetSize()); 
	if (c.GetSize() != a.GetSize()) c.Resize(a.GetSize()); 
#ifdef RV_VECADD
	VectorAdd2_RV(a.GetSize(), a.GetData(), b.GetData(), c.GetData()); 
#else
	for (int i=0; i<a.GetSize(); i++) {
		c[i] = a[i] + b[i]; 
	}
#endif
}

void Subtract(const Vector& a, const Vector& b, Vector& c) {
	CH_TIMERS("vvsub"); 
	CHECK(a.GetSize() == b.GetSize()); 
	if (c.GetSize() != a.GetSize()) c.Resize(a.GetSize()); 
#ifdef RV_VECSUB
	VectorSub2_RV(a.GetSize(), a.GetData(), b.GetData(), c.GetData()); 
#else
	#pragma omp parallel for 
	for (int i=0; i<a.GetSize(); i++) {
		c[i] = a[i] - b[i]; 
	}
#endif
}

void Add(double alpha, const Vector& a, double beta, const Vector& b, Vector& c) {
	CH_TIMERS("vvadd with coeffs"); 
	CHECK(a.GetSize() == b.GetSize()); 
	c.Resize(a.GetSize()); 

	#pragma omp parallel for 
	for (int i=0; i<a.GetSize(); i++) {
		c[i] = alpha*a[i] + beta*b[i]; 
	}
}

} // end namespace fem 