#include "LagrangeSpace.hpp"
#ifdef USE_MPI 
#include <mpi.h>
#endif
#include "Opt.hpp"

using namespace std; 

namespace fem 
{

int FindMatch(const Node& node, const Array<Node>& list) {
	for (int i=0; i<list.GetSize(); i++) {
		if (list[i].CheckEqual(node)) {
			return i; 
		}
	}
	return -1; 
}

int FindMatch(const Node& node, const Array<int>& neighb, 
	Array<Element*> els) {
	for (int e=0; e<neighb.GetSize(); e++) {
		if (neighb[e] >= 0) {
			Element* el = els[neighb[e]]; 
			for (int n=0; n<el->GetNumNodes(); n++) {
				if (el->GetNode(n).CheckEqual(node)) {
					if (el->GetNodeGlobalID(n)>=0)
						return el->GetNodeGlobalID(n); 
				}
			}			
		}
	}
	return -2; 
}

int Search(Element* nel, Node& node, int max) {
	CH_TIMERS("search through neighbors"); 
	for (int i=0; i<nel->GetNumNodes(); i++) {
		if (nel->GetNodeGlobalID(i) >= max) {
			if (nel->GetNode(i).CheckEqual(node)) {
				int id = nel->GetNodeGlobalID(i); 
				if (id>=0) return id; 
			}
		}
	}
	return -1; 
}

bool InVector(const Array<int>& unique, int tag) {
	for (int i=0; i<unique.GetSize(); i++) {
		if (unique[i] == tag) return true; 
	}
	return false; 
}

LagrangeSpace::LagrangeSpace(const Mesh& mesh, int order, int vdim) 
	: FESpace(mesh, order, vdim) {
	CH_TIMERS("setup lagrange space"); 

	// create Elements for every mesh element 
	for (int i=0; i<mesh.GetNumElements(); i++) {
		const MeshEl& el = mesh.GetElement(i); 

		// setup FEM nodes 
		Array<MeshNode> node_list; 
		_mesh.GetNodes(i, node_list); 

		if (el.GetType() == QUAD) {
			_el.Append(new LagrangeQuad(node_list, order)); 						
		} else if (el.GetType() == TRI) {
			_el.Append(new LagrangeTri(node_list, order)); 
		} else if (el.GetType() == HEX) {
			_el.Append(new LagrangeHex(node_list, order)); 
		} else {
			ERROR("element type " << el.GetType() << " not supported"); 
		}

		bool boundary = false; 
		for (int j=0; j<_el[i]->GetNumNodes(); j++) {
			int bc = _el[i]->GetNode(j).GetBC(); 
			// if (bc == DIRICHLET || bc == NEUMANN) {
				// boundary = true; 
			// }
			if (bc != INTERIOR) {
				boundary = true; 
			}
		}
		if (boundary) {
			_bel.Append(_el[i]); 
		}

		if (el.GetNumTags()>0) {
			_el[i]->SetTag(el.GetTag(0)); 
			if (!InVector(_unique_tags, el.GetTag(0))) {
				_unique_tags.Append(el.GetTag(0)); 
			}
		}
	}

	// add linear nodes in 
	{ CH_TIMERS("relabel nodes"); 
	_nodes.Resize(mesh.GetNumNodes()); 
	for (int n=0; n<mesh.GetNumNodes(); n++) {
		const MeshNode& node = mesh.GetNode(n); 
		_nodes[n] = Node(node.x, -1, node.id, node.bc); 
	}

	// add the rest in
	int count = _nodes.GetSize(); 
	for (int e=0; e<_el.GetSize(); e++) {
		Element* el = _el[e]; 
		const MeshEl& mel = mesh.GetElement(e); 
		for (int n=0; n<el->GetNumNodes(); n++) {
			if (el->GetNodeGlobalID(n)<0) {
				Node& node = el->GetNode(n); 
				const Array<int>& nei = mel.GetNeighbors(); 
				int id = -1; 
				for (int ne=0; ne<nei.GetSize(); ne++) {
					if (nei[ne] >= 0 and nei[ne] < e) {
						id = Search(_el[nei[ne]], node, mesh.GetNumNodes()); 
						if (id >= 0) break; 
					}
				}

				if (id>=0) {
					node.SetGlobalID(id); 
				} else {
					node.SetGlobalID(count++); 
					_nodes.Append(node); 
				}
			}
		}
	}
	}

	// assign global id 
	// {
	// 	CH_TIMERS("relabel nodes"); 
	// 	Array<Node> nodes; 
	// 	int count = 0; 
	// 	for (int i=0; i<_el.GetSize(); i++) {
	// 		for (int j=0; j<_el[i]->GetNumNodes(); j++) {
	// 			int match = FindMatch(_el[i]->GetNode(j), nodes); 
	// 			// int match = FindMatch(_el[i]->GetNode(j), 
	// 				// mesh.GetElement(i).GetNeighbors(), _el); 
	// 			if (match >= 0) {
	// 				_el[i]->GetNode(j).SetGlobalID(match); 
	// 			} else {
	// 				_el[i]->GetNode(j).SetGlobalID(count++); 
	// 				nodes.Append(_el[i]->GetNode(j)); 
	// 			}
	// 		}
	// 	}
	// 	_nodes = nodes; 
	// }

	for (int i=0; i<_nodes.GetSize(); i++) {
		if (_nodes[i].GetBC() != INTERIOR) {
			_bnodes.Append(_nodes[i]); 
		}
	}
	_vdofs.Resize(GetNumElements()); 
	for (int e=0; e<GetNumElements(); e++) {
		Element& el = GetEl(e); 
		_vdofs[e].Resize(GetVDim() * el.GetNumNodes()); 
		for (int d=0; d<GetVDim(); d++) {
			for (int i=0; i<el.GetNumNodes(); i++) {
				_vdofs[e][el.GetNumNodes()*d+i] = GetVDim()*el.GetNodeGlobalID(i) + d; 
			}
		}
	}
}

LagrangeLine::LagrangeLine(Array<MeshNode> node_list, int order, int mdim) 
	: Element(node_list, order, LINE, mdim) {

	// store the nodes 
	for (int n=0; n<_geo_nodes.GetSize(); n++) {
		_nodes.Append(Node(_geo_nodes[n].x, _nodes.GetSize(), 
			_geo_nodes[n].id, _geo_nodes[n].bc)); 
		_nodes[n].SetProcs(node_list[n].procs); 
	}

	GenLagrangePolynomials(_order, -1, 1, _p); 
	for (int i=0; i<_p.GetSize(); i++) {
		_dp.Append(_p[i].Derivative()); 
	}
}

void LagrangeLine::CalcShape(Point x, Vector& shape) const {
	shape.SetSize(GetNumNodes()); 
	for (int i=0; i<_p.GetSize(); i++) {
		shape[i] = _p[i](x[0]); 
	}
}

void LagrangeLine::CalcGradShape(Point x, Matrix& gradshape) const {
	gradshape.SetSize(_mdim, GetNumNodes()); 
	for (int d=0; d<_mdim; d++) {
		for (int i=0; i<GetNumNodes(); i++) {
			gradshape(d,i) = _dp[i](x[d]); 
		}
	}
}

LagrangeQuad::LagrangeQuad(Array<MeshNode> node_list, int order, int mdim) 
	: Element(node_list, order, QUAD, mdim) {
	CH_TIMERS("build lagrange quad element"); 

	for (int i=0; i<_geo_nodes.GetSize(); i++) {
		_nodes.Append(Node(_geo_nodes[i].x, _nodes.GetSize(), 
			_geo_nodes[i].id, _geo_nodes[i].bc)); 
		// _nodes.Append(Node(_geo_nodes[i].x, _nodes.GetSize(), 
			// -1, _geo_nodes[i].bc)); 
		_nodes[i].SetProcs(node_list[i].procs); 
	}

	if (_order == 2) {
		CH_TIMERS("setup quadratic"); 
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
			Array<int> procs; 
			_nodes[i].GetProcs().Intersection(_nodes[next].GetProcs(), procs); 
			procs.Sort(); 
			nodes[nodes.GetSize()-1].SetProcs(procs); 
		}
		// add center point 
		Point x; 
		for (int d=0; d<DIM; d++) {
			for (int i=0; i<_nodes.GetSize(); i++) {
				x[d] += _nodes[i].GetX()[d]; 
			}
			x[d] /= _nodes.GetSize(); 
		}
		nodes.Append(Node(x, nodes.GetSize(), -1, INTERIOR)); 
		_nodes = nodes; 
	} 

	else if (_order == 3) {
		CH_TIMERS("setup cubic"); 
		Array<Node> nodes = _nodes; 
		for (int i=0; i<_nodes.GetSize(); i++) {
			int next = (i+1)%_nodes.GetSize(); 
			Point x1; 
			Point x2; 
			for (int d=0; d<DIM; d++) {
				x1[d] = 1./3*(_nodes[next].GetX()[d] - 
					_nodes[i].GetX()[d]) + _nodes[i].GetX()[d]; 
				x2[d] = 2./3*(_nodes[next].GetX()[d] - 
					_nodes[i].GetX()[d]) + _nodes[i].GetX()[d]; 
			}

			int bc; 
			if (_nodes[i].GetBC() == nodes[next].GetBC()) {
				bc = _nodes[i].GetBC(); 
			} else {
				bc = INTERIOR; 
			}
			Array<int> procs; 
			_nodes[i].GetProcs().Intersection(_nodes[next].GetProcs(), procs); 
			procs.Sort(); 
			nodes.Append(Node(x1, nodes.GetSize(), -1, bc)); 
			nodes.Append(Node(x2, nodes.GetSize(), -1, bc)); 
			nodes[nodes.GetSize()-1].SetProcs(procs); 
			nodes[nodes.GetSize()-2].SetProcs(procs); 
		}

		nodes.Append(Node(
			{nodes[4].GetX()[0], nodes[11].GetX()[1]}, 
			nodes.GetSize(), -1, INTERIOR)); 
		nodes.Append(Node(
			{nodes[5].GetX()[0], nodes[11].GetX()[1]}, 
			nodes.GetSize(), -1, INTERIOR)); 
		nodes.Append(Node(
			{nodes[4].GetX()[0], nodes[10].GetX()[1]}, 
			nodes.GetSize(), -1, INTERIOR)); 
		nodes.Append(Node(
			{nodes[5].GetX()[0], nodes[10].GetX()[1]}, 
			nodes.GetSize(), -1, INTERIOR)); 
		_nodes = nodes; 
	}

	GenLagrangePolynomials(_order, -1, 1, _p); 
	if (_order == 1) {
		CH_TIMERS("setup linear shape"); 
		_pp.Append(PolyProduct(_p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[1])); 
		_pp.Append(PolyProduct(_p[0], _p[1])); 
	} else if (_order == 2) {
		CH_TIMERS("setup quadratic shape"); 
		_pp.Append(PolyProduct(_p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[2], _p[0]));
		_pp.Append(PolyProduct(_p[2], _p[2]));
		_pp.Append(PolyProduct(_p[0], _p[2]));

		_pp.Append(PolyProduct(_p[1], _p[0]));
		_pp.Append(PolyProduct(_p[2], _p[1]));
		_pp.Append(PolyProduct(_p[1], _p[2]));
		_pp.Append(PolyProduct(_p[0], _p[1]));
		_pp.Append(PolyProduct(_p[1], _p[1]));
	} else if (_order == 3) {
		CH_TIMERS("setup cubic shape"); 
		// 0-3
		_pp.Append(PolyProduct(_p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[3], _p[0])); 
		_pp.Append(PolyProduct(_p[3], _p[3])); 
		_pp.Append(PolyProduct(_p[0], _p[3])); 

		// 4 - 11 
		_pp.Append(PolyProduct(_p[1], _p[0])); // 4 
		_pp.Append(PolyProduct(_p[2], _p[0])); // 5 
		_pp.Append(PolyProduct(_p[3], _p[1])); // 6 
		_pp.Append(PolyProduct(_p[3], _p[2])); // 7 
		_pp.Append(PolyProduct(_p[2], _p[3])); // 8 
		_pp.Append(PolyProduct(_p[1], _p[3])); // 9 
		_pp.Append(PolyProduct(_p[0], _p[2])); // 10 
		_pp.Append(PolyProduct(_p[0], _p[1])); // 11 

		// 12 - 15
		_pp.Append(PolyProduct(_p[1], _p[1])); // 12
		_pp.Append(PolyProduct(_p[2], _p[1])); // 13 
		_pp.Append(PolyProduct(_p[1], _p[2])); // 14 
		_pp.Append(PolyProduct(_p[2], _p[2])); // 15
	} else {
		ERROR("order " << _order << " not defined"); 
	}

	{CH_TIMERS("gradients"); 
	// for (int i=0; i<_p.GetSize(); i++) {
	// 	_dp.Append(_p[i].Derivative()); 
	// }

	_dpp.Resize(_pp.GetSize()); 
	for (int i=0; i<_pp.GetSize(); i++) {
		_pp[i].Gradient(_dpp[i]); 
	}
	}

	{CH_TIMERS("setup batched"); 
	_shapex.Resize(_pp.GetSize()*(_order+1)); 
	_shapey.Resize(_pp.GetSize()*(_order+1)); 

	for (int i=0; i<_pp.GetSize(); i++) {
		for (int j=0; j<_order+1; j++) {
			_shapex[i*(_order+1)+j] = _pp[i].Get1D(0).GetCoefs()[j];
			_shapey[i*(_order+1)+j] = _pp[i].Get1D(1).GetCoefs()[j];  
		}
	}

	_dshapex.Resize(_mdim*GetNumNodes()*(_order+1));
	_dshapey.Resize(_mdim*GetNumNodes()*(_order+1));  
	for (int i=0; i<GetNumNodes(); i++) {
		for (int j=0; j<_order; j++) {
			_dshapex[i*(_order+1)+j] = _dpp[i][0].Get1D(0).GetCoefs()[j]; 
		} 
		for (int j=0; j<_order+1; j++) {
			_dshapey[i*(_order+1)+j] = _dpp[i][0].Get1D(1).GetCoefs()[j]; 
		}
	}
	for (int i=0; i<GetNumNodes(); i++) {
		for (int j=0; j<_order+1; j++) {
			_dshapex[i*(_order+1)+j+GetNumNodes()*(_order+1)] = 
				_dpp[i][1].Get1D(0).GetCoefs()[j]; 
		} 
		for (int j=0; j<_order; j++) {
			_dshapey[i*(_order+1)+j+GetNumNodes()*(_order+1)] = 
				_dpp[i][1].Get1D(1).GetCoefs()[j]; 
		}
	}
	}
}

void LagrangeQuad::CalcShape(Point x, 
	Vector& shape) const {
	CH_TIMERS("lagrange quad calc shape"); 
	shape.SetSize(GetNumNodes()); 
#ifdef RV_SHAPE 
	shape = 1.; 
	CalcShape_RV(GetNumNodes(), _order+1, _shapex.GetData(), 
		_shapey.GetData(), shape.GetData(), x.GetData()); 
#else
	for (int i=0; i<_pp.GetSize(); i++) {
		shape[i] = _pp[i](x); 
	}
#endif
}

void LagrangeQuad::CalcGradShape(Point x, 
	Matrix& gradshape) const {
	CH_TIMERS("lagrange quad calc grad shape"); 
	gradshape.SetSize(_mdim, GetNumNodes()); 
#ifdef RV_GSHAPE 
	gradshape = 1.; 
	CalcShape_RV(GetNumNodes()*2, _order+1, _dshapex.GetData(), 
		_dshapey.GetData(), gradshape.GetData(), x.GetData()); 
#else
	for (int i=0; i<_dpp.GetSize(); i++) {
		for (int j=0; j<_mdim; j++) {
			gradshape(j, i) = _dpp[i][j](x); 
		}
	}
#endif
}

void LagrangeQuad::CalcPhysGradShape(ElTrans& trans, 
	Matrix& pgshape) const {
	Matrix gshape; 
	CalcGradShape(trans.GetX(), gshape); 
	trans.InverseJacobian().Mult(gshape, pgshape); 
}

LagrangeTri::LagrangeTri(Array<MeshNode> node_list, int order, int mdim) 
	: Element(node_list, order, TRI, mdim) {

	for (int i=0; i<_geo_nodes.GetSize(); i++) {
		_nodes.Append(Node(_geo_nodes[i].x, _nodes.GetSize(), 
			_geo_nodes[i].id, _geo_nodes[i].bc)); 
		_nodes[i].SetProcs(_geo_nodes[i].procs); 
	}

	if (_order == 1) {
		_p.Append(Poly1D({1,-1})); 
		_p.Append(Poly1D({0,1})); 
	}

	else {
		ERROR("order " << _order << " not defined"); 
	}

	for (int i=0; i<_p.GetSize(); i++) {
		_dp.Append(_p[i].Derivative()); 
	}
}

void LagrangeTri::CalcShape(Point x, Vector& shape) const {
	shape.SetSize(GetNumNodes()); 

	shape[0] = 1. - x[0] - x[1]; 
	shape[1] = x[0]; 
	shape[2] = x[1]; 
}

void LagrangeTri::CalcGradShape(Point x, Matrix& gradshape) const {
	gradshape.SetSize(_dim, GetNumNodes()); 

	gradshape(0,0) = -1.;  
	gradshape(1,0) = -1.; 

	gradshape(0,1) = 1.; 
	gradshape(1,1) = 0.; 

	gradshape(0,2) = 0.; 
	gradshape(1,2) = 1.; 
}

void LagrangeTri::CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const {
	Matrix gshape; 
	CalcGradShape(trans.GetX(), gshape); 
	trans.InverseJacobian().Mult(gshape, pgshape); 
}

LagrangeHex::LagrangeHex(Array<MeshNode> node_list, int order, int mdim) : Element(node_list, order, HEX, mdim) {

	// append nodes 
	for (int i=0; i<_geo_nodes.GetSize(); i++) {
		_nodes.Append(Node(_geo_nodes[i].x, _nodes.GetSize(), _geo_nodes[i].id, _geo_nodes[i].bc)); 
		_nodes[i].SetProcs(node_list[i].procs); 
	}

	if (_order == 2) {
		Array<Node> nodes = _nodes; 

		// bottom layer 
		for (int i=0; i<4; i++) {
			int next = (i+1)%4; 
			Point x; 
			for (int d=0; d<DIM; d++) {
				x[d] = .5*(_nodes[i].GetX()[d] + _nodes[next].GetX()[d]); 
			}

			int bc; 
			if (_nodes[i].GetBC()==_nodes[next].GetBC()) {
				bc = _nodes[i].GetBC(); 
			} else bc = INTERIOR; 
			nodes.Append(Node(x, _nodes.GetSize()+i, -1, bc)); 
			Array<int> procs; 
			_nodes[i].GetProcs().Intersection(_nodes[next].GetProcs(), procs); 
			procs.Sort(); 
			nodes[nodes.GetSize()-1].SetProcs(procs); 
		}

		// add the center point 
		Point x; 
		for (int d=0; d<DIM; d++) {
			for (int i=0; i<4; i++) {
				x[d] += _nodes[i].GetX()[d]; 				
			}
			x[d] /= 4; 
		}
		nodes.Append(Node(x, nodes.GetSize(), -1, INTERIOR)); 

		// middle layer 
		Array<Node> corners; 
		for (int i=0; i<4; i++) {
			int next = i+4; 
			x = {0,0,0}; 
			for (int d=0; d<DIM; d++) {
				x[d] = .5*(_nodes[i].GetX()[d] + _nodes[next].GetX()[d]); 
			}
			int bc; 
			if (_nodes[i].GetBC() == _nodes[next].GetBC()) bc = _nodes[i].GetBC(); 
			else bc = INTERIOR; 
			nodes.Append(Node(x, nodes.GetSize(), -1, bc)); 
			corners.Append(nodes[nodes.GetSize()-1]); 
		}

		for (int i=0; i<4; i++) {
			int next = (i+1)%4; 
			x = {0,0,0}; 
			for (int d=0; d<DIM; d++) {
				x[d] = .5*(corners[i].GetX()[d] + corners[next].GetX()[d]); 
			}
			int bc; 
			if (corners[i].GetBC() == corners[next].GetBC()) bc = corners[i].GetBC(); 
			else bc = INTERIOR; 
			nodes.Append(Node(x, nodes.GetSize(), -1, bc)); 
		}
		// add center point 
		x = {0,0,0}; 
		for (int d=0; d<DIM; d++) {
			for (int i=0; i<4; i++) {
				x[d] += corners[i].GetX()[d]; 
			}
			x[d] /= 4; 
		}
		nodes.Append(Node(x, nodes.GetSize(), -1, INTERIOR)); 

		// top layer 
		for (int i=4; i<8; i++) {
			int next = (i+1)%4 + 4; 
			Point x; 
			for (int d=0; d<DIM; d++) {
				x[d] = .5*(_nodes[i].GetX()[d] + _nodes[next].GetX()[d]); 
			}

			int bc; 
			if (_nodes[i].GetBC()==_nodes[next].GetBC()) {
				bc = _nodes[i].GetBC(); 
			} else bc = INTERIOR; 
			nodes.Append(Node(x, nodes.GetSize(), -1, bc)); 
			Array<int> procs; 
			_nodes[i].GetProcs().Intersection(_nodes[next].GetProcs(), procs); 
			procs.Sort(); 
			nodes[nodes.GetSize()-1].SetProcs(procs); 
		}
		// add center point 
		x = {0,0,0}; 
		for (int d=0; d<DIM; d++) {
			for (int i=4; i<8; i++) {
				x[d] += _nodes[i].GetX()[d]; 
			}
			x[d] /= 4; 
		}
		nodes.Append(Node(x, nodes.GetSize(), -1, INTERIOR)); 

		_nodes = nodes; 
	} else if (_order > 2) ERROR("order " << _order << " not implemented"); 

	GenLagrangePolynomials(_order, -1, 1, _p); 
	if (_order == 1) {
		// bottom quad 
		_pp.Append(PolyProduct(_p[0], _p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[1], _p[1], _p[0])); 
		_pp.Append(PolyProduct(_p[0], _p[1], _p[0])); 

		// top quad 
		_pp.Append(PolyProduct(_p[0], _p[0], _p[1])); 
		_pp.Append(PolyProduct(_p[1], _p[0], _p[1])); 
		_pp.Append(PolyProduct(_p[1], _p[1], _p[1])); 
		_pp.Append(PolyProduct(_p[0], _p[1], _p[1])); 
	} else if (_order == 2) {
		// bottom four 
		_pp.Append(PolyProduct(_p[0], _p[0], _p[0])); 
		_pp.Append(PolyProduct(_p[2], _p[0], _p[0]));
		_pp.Append(PolyProduct(_p[2], _p[2], _p[0]));
		_pp.Append(PolyProduct(_p[0], _p[2], _p[0]));
		// top four 
		_pp.Append(PolyProduct(_p[0], _p[0], _p[2])); 
		_pp.Append(PolyProduct(_p[2], _p[0], _p[2]));
		_pp.Append(PolyProduct(_p[2], _p[2], _p[2]));
		_pp.Append(PolyProduct(_p[0], _p[2], _p[2]));
		
		// bottom ho nodes 
		_pp.Append(PolyProduct(_p[1], _p[0], _p[0]));
		_pp.Append(PolyProduct(_p[2], _p[1], _p[0]));
		_pp.Append(PolyProduct(_p[1], _p[2], _p[0]));
		_pp.Append(PolyProduct(_p[0], _p[1], _p[0]));
		_pp.Append(PolyProduct(_p[1], _p[1], _p[0]));

		// middle ho nodes
		_pp.Append(PolyProduct(_p[0], _p[0], _p[1])); 
		_pp.Append(PolyProduct(_p[2], _p[0], _p[1]));
		_pp.Append(PolyProduct(_p[2], _p[2], _p[1]));
		_pp.Append(PolyProduct(_p[0], _p[2], _p[1]));

		_pp.Append(PolyProduct(_p[1], _p[0], _p[1]));
		_pp.Append(PolyProduct(_p[2], _p[1], _p[1]));
		_pp.Append(PolyProduct(_p[1], _p[2], _p[1]));
		_pp.Append(PolyProduct(_p[0], _p[1], _p[1]));
		_pp.Append(PolyProduct(_p[1], _p[1], _p[1]));
		
		// top ho nodes 
		_pp.Append(PolyProduct(_p[1], _p[0], _p[2]));
		_pp.Append(PolyProduct(_p[2], _p[1], _p[2]));
		_pp.Append(PolyProduct(_p[1], _p[2], _p[2]));
		_pp.Append(PolyProduct(_p[0], _p[1], _p[2]));
		_pp.Append(PolyProduct(_p[1], _p[1], _p[2]));
	}
	CHECK(_pp.GetSize() == _nodes.GetSize()); 

	// get gradients 
	_dpp.Resize(_pp.GetSize()); 
	for (int i=0; i<_pp.GetSize(); i++) {
		_pp[i].Gradient(_dpp[i]); 
	}
}

void LagrangeHex::CalcShape(Point x, Vector& shape) const {
	shape.Resize(GetNumNodes()); 
	for (int i=0; i<GetNumNodes(); i++) {
		shape[i] = _pp[i](x); 
	}
}

void LagrangeHex::CalcGradShape(Point x, Matrix& gradshape) const {
	gradshape.SetSize(_mdim, GetNumNodes()); 
	for (int i=0; i<_dpp.GetSize(); i++) {
		for (int j=0; j<_mdim; j++) {
			gradshape(j, i) = _dpp[i][j](x); 
		}
	}
}

void LagrangeHex::CalcPhysGradShape(ElTrans& trans, Matrix& pgshape) const {
	Matrix gshape; 
	CalcGradShape(trans.GetX(), gshape); 
	trans.InverseJacobian().Mult(gshape, pgshape); 
}

} // end namespace fem 