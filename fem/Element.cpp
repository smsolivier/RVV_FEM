#include "Element.hpp" 

using namespace std; 

namespace fem 
{

Element::Element(Array<MeshNode>& node_list, int order, int geo, int mdim) {
	_geo_nodes = node_list; 
	_order = order; 
	_dim = GetGeoDim(geo); 
	_mdim = mdim; 
	_geo = geo; 
	_volume = -1; 

	_neighbors.Resize(FacesPerEl(GetType())); 
	_neighbors = -1; 

	_faceRefTrans.Resize(_neighbors.GetSize()); 
	_faceRefTrans = NULL; 

	_faceTrans.Resize(_neighbors.GetSize()); 
	_faceTrans = NULL; 

	_trans = new ElTrans(*this); 
}

Element::~Element() {
	if (_trans) delete _trans; 
}

Node& Element::GetNode(int index) {
	return _nodes[index]; 
}

const Node& Element::GetNode(int index) const {
	return _nodes[index]; 
}

int Element::GetNodeGlobalID(int i) const {
	return _nodes[i].GetGlobalID(); 
}

const Matrix& Element::GetNodeLocationMatrix() {
	if (_nodeLocMatrix.Height() == 0) {
		_nodeLocMatrix.SetSize(GetNumNodes(), GetDim()); 
		for (int n=0; n<GetNumNodes(); n++) {
			for (int d=0; d<GetDim(); d++) {
				_nodeLocMatrix(n,d) = GetNode(n).GetX()[d]; 
			}
		}
	}
	return _nodeLocMatrix; 
}

ElTrans& Element::GetFaceRefTrans(int e) {
	int ind = FindNeighbor(e); 
	CHECKMSG(ind >= 0, "neighbor not found"); 
	CHECKMSG(_faceRefTrans[ind], "reference face transformation not initialized"); 
	return _faceRefTrans[ind]->GetTrans(); 
}

ElTrans& Element::GetFaceRefTransFromFace(int f) {
	return _faceRefTrans[f]->GetTrans(); 
}

ElTrans& Element::GetFaceTrans(int e) {
	int ind = FindNeighbor(e); 
	CHECKMSG(_faceTrans[ind], "face transformation not initialized"); 
	return _faceTrans[ind]->GetTrans(); 
}

ElTrans& Element::GetFaceTransFromFace(int f) {
	return _faceTrans[f]->GetTrans(); 
}

int Element::FindNeighbor(int e) const {
	for (int i=0; i<_neighbors.GetSize(); i++) {
		if (_neighbors[i] == e) {
			return i; 
		}
	}
	return -1; 
}

void Element::CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const {
	Matrix gshape; 
	CalcGradShape(trans.GetX(), gshape); 
	trans.InverseJacobian().Mult(gshape, pgshape); 
}

double Element::Integrate(const Vector& x, Quadrature* quad) const {
	if (quad == NULL) {
		quad = QRules.Get(GetType(), INTEGRATION_ORDER, INTEGRATION_TYPE); 
	}

	double sum = 0.; 
	Vector shape; 
	for (int i=0; i<quad->NumPoints(); i++) {
		_trans->SetX(quad->X(i)); 
		CalcShape(quad->X(i), shape); 
		double eval = 0.; 
		for (int j=0; j<GetNumNodes(); j++) {
			eval += shape[j] * x[GetNodeGlobalID(j)]; 
		}
		sum += eval * quad->Weight(i) * _trans->Determinant(); 
	}
	return sum; 
}

Point Element::Centroid() const {
	Point c; 
	for (int i=0; i<GetNumNodes(); i++) {
		c += GetNode(i).GetX(); 
	}
	c /= GetNumNodes(); 
	Point phys; 
	_trans->Transform(c, phys); 
	return phys; 
}

double Element::Volume(Quadrature* quad) {
	if (_volume != -1) return _volume; 

	if (quad == NULL) {
		quad = QRules.Get(INTEGRATION_ORDER, GetType(), INTEGRATION_TYPE); 
	} 

	_volume = 0.; 
	for (int i=0; i<quad->NumPoints(); i++) {
		_trans->SetX(quad->X(i)); 
		_volume += quad->Weight(i) * _trans->Determinant(); 
	}

	return _volume; 
}

double Element::EnergyNorm(const Vector& x, Quadrature* quad) const {
	if (quad == NULL) {
		quad = QRules.Get(INTEGRATION_ORDER, GetType(), INTEGRATION_TYPE); 
	}

	double e = 0; 
	Matrix gradshape; 
	Vector grad; 
	for (int i=0; i<quad->NumPoints(); i++) {
		_trans->SetX(quad->X(i)); 
		double val = pow(Interpolate(x, quad->X(i)), 2); 
		InterpolateGradient(x, quad->X(i), grad); 
		CHECKMSG(grad.Dot(grad) >= 0, "nabla dot nabla = " << grad.Dot(grad)); 
		e += (grad.Dot(grad) + val) * quad->Weight(i) * _trans->Determinant(); 
		// e += val * quad->Weight(i) * trans.Determinant(); 
	}
	CHECK(isfinite(e)); 
	return e; 
}

double Element::Interpolate(const Vector& x, const Point& p) const {
	CHECK(x.GetSize() > 0); 
	double sum = 0; 
	Vector shape; 
	CalcShape(p, shape); 
	for (int i=0; i<GetNumNodes(); i++) {
		sum += x[GetNodeGlobalID(i)] * shape[i]; 
	}
	return sum; 
}

void Element::InterpolateGradient(const Vector& x, 
	const Point& p, Vector& grad) const {

	CHECK(x.GetSize()>0); 
	_trans->SetX(p); 
	grad.SetSize(_dim); 
	grad = 0.; 
	Matrix gradshape; 
	CalcPhysGradShape(*_trans, gradshape); 
	for (int i=0; i<_dim; i++) {
		for (int j=0; j<GetNumNodes(); j++) {
			grad[i] += gradshape(i,j) * x[GetNodeGlobalID(j)]; 
		}
	}
}

void Element::Print(std::ostream& out) const {
	out << "Element Information: " << endl; 
	out << "\tElement Geometry = " << GetType() << endl; 
	out << "\tNumber of Nodes = " << GetNumNodes() << endl; 
	out << "\tFEM order = " << _order << endl; 
	out << "\tNode IDs = "; 
	for (int i=0; i<GetNumNodes(); i++) {
		out << _nodes[i].GetGlobalID() << " "; 
	}
	out << endl; 
	out << "\tNeighbors = "; 
	for (int i=0; i<_neighbors.GetSize(); i++) {
		out << _neighbors[i] << " "; 
	}
	out << endl; 
}

} // end namespace fem 