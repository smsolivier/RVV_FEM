#include "L2Space.hpp"

using namespace std; 

namespace fem 
{

L2Space::L2Space(const Mesh& mesh, int order) : FESpace(mesh, order) {

	// assign FEM nodes by mesh elements 
	for (int i=0; i<mesh.GetNumElements(); i++) {
		const MeshEl& el = mesh.GetElement(i); 
		ASSERT(el.GetType()==QUAD); 

		Array<MeshNode> node_list(4); 
		for (int j=0; j<4; j++) {
			node_list[j] = mesh.GetNode(el[j]); 
		}
		_el.Append(new L2Quad(node_list, _order)); 
		_el[i]->SetID(i); 

		const Array<int>& neighb = el.GetNeighbors(); 
		for (int j=0; j<neighb.GetSize(); j++) {
			_el[i]->SetNeighbor(j, neighb[j]); 
		}
	}

	// set global id's of elements 
	int count = 0; 
	for (int i=0; i<_el.GetSize(); i++) {
		for (int j=0; j<_el[i]->GetNumNodes(); j++) {
			_el[i]->GetNode(j).SetGlobalID(count++); 
		}
	}

	// store nodes 
	for (int i=0; i<_el.GetSize(); i++) {
		for (int j=0; j<_el[i]->GetNumNodes(); j++) {
			_nodes.Append(_el[i]->GetNode(j)); 
		}
	}
}

L2Segment::L2Segment(Array<MeshNode> list, int order, int mdim) 
	: Element(list, order, LINE, mdim) {
	if (_order == 0) {
		Point avg; 
		for (int i=0; i<list.GetSize(); i++) {
			avg += list[i].x; 
		}
		avg /= list.GetSize(); 
		_nodes.Append(Node(avg, 0, -1, INTERIOR)); 
	} else if (_order == 1) {
		for (int i=0; i<list.GetSize(); i++) {
			_nodes.Append(Node(list[i].x, i, -1, list[i].bc)); 
		}
	} else if (_order == 2) {
		for (int i=0; i<list.GetSize(); i++) {
			_nodes.Append(Node(list[i].x, i, -1, list[i].bc)); 
		}
		Point mid; 
		mid += list[0].x; 
		mid += list[1].x; 
		mid /= 2; 
		int bc = (list[0].bc == list[1].bc) ? list[0].bc : INTERIOR; 
		_nodes.Append(Node(mid, 2, -1, bc)); 
	} else {
		ERROR("order " << _order << " not implemented for L2Segment"); 
	}

	GenLagrangePolynomials(_order, -1, 1, _p); 
	// if (_order == 0) {
	// 	_p.Append(Poly1D({1})); 
	// } else if (_order == 1) {
	// 	_p.Append(Poly1D({1, -1})); 
	// 	_p.Append(Poly1D({0, 1})); 
	// } else if (_order == 2) {
	// 	_p.Append(Poly1D({1,-3,2})); 
	// 	_p.Append(Poly1D({0,-1,2})); 
	// 	_p.Append(Poly1D({0,4,-4}));
	// }

	for (int i=0; i<_p.GetSize(); i++) {
		_dp.Append(_p[i].Derivative()); 
	}
}

void L2Segment::CalcShape(Point x, Vector& shape) const {
	shape.SetSize(GetNumNodes()); 
	for (int i=0; i<_p.GetSize(); i++) {
		shape[i] = _p[i](x[0]); 
	}
}

void L2Segment::CalcGradShape(Point x, Matrix& gradshape) const {
	gradshape.SetSize(_mdim, GetNumNodes()); 
	for (int d=0; d<_mdim; d++) {
		for (int i=0; i<GetNumNodes(); i++) {
			gradshape(d,i) = _dp[i](x[d]); 
		}
	}
}

L2Quad::L2Quad(Array<MeshNode> node_list, int order, int mdim) 
	: Element(node_list, order, QUAD, mdim) {
	if (_order == 0) {
		Point avg; 
		for (int i=0; i<node_list.GetSize(); i++) {
			avg += node_list[i].x; 
		}
		avg /= node_list.GetSize(); 
		_nodes.Append(Node(avg, 0, -1, INTERIOR));
	}
	else if (_order > 0) {
		for (int i=0; i<node_list.GetSize(); i++) {
			_nodes.Append(Node(node_list[i].x, i, -1, node_list[i].bc));  
		}
	} 

	if (_order == 0) {

	} 
	else if (_order == 1) {

	}
	else if (_order == 2) {
		Array<Node> nodes = _nodes; 
		for (int i=0; i<_nodes.GetSize(); i++) {
			int next = (i+1)%_nodes.GetSize(); 
			Point x; 
			for (int d=0; d<DIM; d++) {
				x[d] = .5*(_nodes[i].GetX()[d] + _nodes[next].GetX()[d]); 
			}

			int bc; 
			if (_nodes[i].GetBC() == _nodes[next].GetBC()) {
				bc = _nodes[i].GetBC(); 
			} else {
				bc = INTERIOR; 
			}
			nodes.Append(Node(x, _nodes.GetSize()+i, -1, bc)); 
		}
		// add center point 
		Point x; 
		for (int d=0; d<DIM; d++) {
			for (int i=0; i<_nodes.GetSize(); i++) {
				x[d] += _nodes[i].GetX()[d]; 
			}
		}
		x /= _nodes.GetSize(); 
		nodes.Append(Node(x, nodes.GetSize(), -1, INTERIOR)); 
		_nodes = nodes; 
	}
	else {
		ERROR("order " << _order << " not implemented for L2"); 
	}

	if (_order == 0) {
		_p.Append(Poly1D({1})); 
		_pp.Append(PolyProduct(_p[0])); 
	}
	else if (_order == 1) {
		_p.Append(Poly1D({1, -1})); 
		_p.Append(Poly1D({0, 1})); 
		_pp.Append(PolyProduct(_p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[1])); 
		_pp.Append(PolyProduct(_p[0], _p[1])); 
	} else if (_order == 2) {
		_p.Append(Poly1D({1,-3,2})); 
		_p.Append(Poly1D({0,-1,2})); 
		_p.Append(Poly1D({0,4,-4}));

		_pp.Append(PolyProduct(_p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[0]));
		_pp.Append(PolyProduct(_p[1], _p[1]));
		_pp.Append(PolyProduct(_p[0], _p[1]));

		_pp.Append(PolyProduct(_p[2], _p[0]));
		_pp.Append(PolyProduct(_p[1], _p[2]));
		_pp.Append(PolyProduct(_p[2], _p[1]));
		_pp.Append(PolyProduct(_p[0], _p[2]));
		_pp.Append(PolyProduct(_p[2], _p[2]));
	}

	// get derivatives of _p 
	for (int i=0; i<_p.GetSize(); i++) {
		_dp.Append(_p[i].Derivative()); 
	}

	_dpp.Resize(_pp.GetSize()); 
	for (int i=0; i<_pp.GetSize(); i++) {
		_pp[i].Gradient(_dpp[i]); 		
	}

	// setup face transformations 
	{
		MeshNode n0 = {0, {0,0}, INTERIOR}; 
		MeshNode n1 = {0, {1,0}, INTERIOR}; 
		_faceRefTrans[0] = new L2Segment({n0, n1}, _order, _mdim); 
	}
	{
		MeshNode n0 = {0, {1,0}, INTERIOR}; 
		MeshNode n1 = {0, {1,1}, INTERIOR}; 
		_faceRefTrans[1] = new L2Segment({n0, n1}, _order, _mdim); 
	}
	{
		MeshNode n0 = {0, {1,1}, INTERIOR}; 
		MeshNode n1 = {0, {0,1}, INTERIOR}; 
		_faceRefTrans[2] = new L2Segment({n0, n1}, _order, _mdim); 
	}
	{
		MeshNode n0 = {0, {0,1}, INTERIOR}; 
		MeshNode n1 = {0, {0,0}, INTERIOR}; 
		_faceRefTrans[3] = new L2Segment({n0, n1}, _order, _mdim); 
	}
	// setup reference to physical space face transformations 
	for (int i=0; i<_faceTrans.GetSize(); i++) {
		int next = (i+1)%_faceTrans.GetSize(); 
		_faceTrans[i] = new L2Segment({_geo_nodes[i], _geo_nodes[next]}, 
			_order, _mdim); 
	}
}

void L2Quad::CalcShape(Point x, Vector& shape) const {
	shape.SetSize(GetNumNodes()); 

	for (int i=0; i<_pp.GetSize(); i++) {
		shape[i] = _pp[i](x); 
	}
}

void L2Quad::CalcGradShape(Point x, Matrix& gradshape) const {
	gradshape.SetSize(_mdim, GetNumNodes()); 

	for (int i=0; i<_dpp.GetSize(); i++) {
		for (int j=0; j<_mdim; j++) {
			gradshape(j, i) = _dpp[i][j](x); 
		}
	}
}

} // end namespace fem 