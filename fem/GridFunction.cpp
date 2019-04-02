#include "GridFunction.hpp"
#ifdef USE_MPI 
#include <mpi.h>
#endif

namespace fem 
{

GridFunction::GridFunction() : Vector() {
	_space = NULL; 
}

GridFunction::GridFunction(const FESpace* space) : Vector(space->GetVSize(), 0) {
	_space = space; 
}

void GridFunction::Project(double (*f)(const Point&)) {
	// for (int i=0; i<_space->GetNumNodes(); i++) {
	// 	(*this)[i] = f(_space->GetNode(i).GetX()); 
	// }

	for (int i=0; i<_space->GetNumElements(); i++) {
		Array<int> vdofs; 
		_space->GetVDofs(i, vdofs); 
		Element& el = _space->GetEl(i); 
		for (int j=0; j<vdofs.GetSize(); j++) {
			(*this)[vdofs[j]] = f(el.GetNode(j).GetX()); 
		}
	}
}

double GridFunction::EnergyNorm(Quadrature* quad) const {
	double e = 0; 
	for (int i=0; i<_space->GetNumElements(); i++) {
		const Element& el = _space->GetEl(i); 
		quad = QRules.Get(el.GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
		e += el.EnergyNorm(*this, quad); 
	}

	return e; 
}

double GridFunction::L2Error(Coefficient* c) const {
	double val; 
	double dot; 
	Vector shape; 
	Array<int> vdofs; 
	Vector vals; 
	double error = 0.; 
	for (int n=0; n<_space->GetNumElements(); n++) {
		Element& el = _space->GetEl(n); 
		Quadrature* quad = QRules.Get(el.GetType(), 
			INTEGRATION_ORDER, INTEGRATION_TYPE); 
		ElTrans& trans = el.GetTrans(); 
		_space->GetVDofs(n, vdofs); 
		for (int i=0; i<quad->NumPoints(); i++) {
			Point ip = quad->X(i); 
			val = c->Eval(trans, ip); 
			el.CalcShape(ip, shape); 
			GetFromDofs(vdofs, vals); 
			dot = vals * shape; 
			error += quad->Weight(i) * pow(dot - val, 2) * fabs(trans.Weight()); 
		}
	}
	return sqrt(error); 
}

} // end namespace fem 