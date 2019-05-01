#include "Quadrature.hpp"
#include "Mesh.hpp"

#define QUAD_ARRAY_SIZE 20

using namespace std; 

namespace fem 
{

QuadRules QRules; 

QuadLegendre::QuadLegendre(int p) : Quadrature(p) {
	Array<double> x; 
	switch(_p) {
		case 1: {
			x = {.5}; 
			_w = {1}; 
			break; 
		}
		case 2: {
			x = {0.21132487, 0.78867513}; 
			_w = {.5, .5}; 
			break; 
		}

		case 3 : {
			x = {0.11270167, 0.5, 0.88729833}; 
			_w = {0.27777778, 0.44444444, 0.27777778}; 
			break; 
		}

		case 4: {
			x = {0.06943184, 0.33000948, 0.66999052, 0.93056816}; 
			_w = {0.17392742, 0.32607258, 0.32607258, 0.17392742}; 
			break; 
		}

		case 5: {
			x = {0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992}; 
			_w = {0.11846344, 0.23931434, 0.28444444, 0.23931434, 0.11846344}; 
			break; 
		}

		case 6: {
			x = {0.03376524, 0.16939531, 0.38069041, 0.61930959, 
				0.83060469, 0.96623476}; 
			_w = {0.08566225, 0.18038079, 0.23395697, 0.23395697, 
				0.18038079, 0.08566225}; 
			break; 
		}

		case 7: {
			x = {0.02544604, 0.12923441, 0.29707742, 0.5, 
				0.70292258, 0.87076559, 0.97455396}; 
			_w = {0.06474248, 0.1398527 , 0.19091503, 0.20897959, 
				0.19091503, 0.1398527 , 0.06474248}; 
			break; 
		}

		case 8: {
			x = {0.01985507, 0.10166676, 0.2372338 , 0.40828268, 0.59171732, 
				0.7627662 , 0.89833324, 0.98014493}; 
			_w = {0.05061427, 0.11119052, 0.15685332, 0.18134189, 0.18134189, 
				0.15685332, 0.11119052, 0.05061427}; 
			break; 
		}

		case 9: {
			x = {0.01591988, 0.08198445, 0.19331428, 0.33787329, 0.5,
				0.66212671, 0.80668572, 0.91801555, 0.98408012}; 
			_w = {0.04063719, 0.09032408, 0.13030535, 0.15617354, 0.16511968, 
				0.15617354, 0.13030535, 0.09032408, 0.04063719}; 
			break; 
		}

		case 10: {
			x = {0.01304674, 0.06746832, 0.16029522, 0.2833023 , 0.42556283, 
				0.57443717, 0.7166977 , 0.83970478, 0.93253168, 0.98695326}; 
			_w = {0.03333567, 0.07472567, 0.10954318, 0.13463336, 0.14776211, 
				0.14776211, 0.13463336, 0.10954318, 0.07472567, 0.03333567}; 
			break; 
		}

		default : {
			ERROR("legendre integration of order " << _p 
				<< " not implemented"); 
		}
	}

	for (int i=0; i<x.GetSize(); i++) {
		_x.Append(Point({(x[i] - .5)*2.})); 
		_w[i] *= 2.; 
	}
}

QuadLobatto::QuadLobatto(int p) : Quadrature(p) {
	Array<double> x; 
	switch(_p) {
		case 2: {
			x = {-1., 1.};
			_w = {1, 1};  
			break; 
		}

		case 3: {
			x = {-1., 0., 1.}; 
			_w = {1./3, 4./3, 1./3}; 
			break; 
		}

		case 4: {
			x = {-1., -1./5*sqrt(5), 1./5*sqrt(5), 1.}; 
			_w = {1./6, 5./6, 5./6, 1./6}; 
			break; 
		}

		case 5: {
			x = {0, 1./7*sqrt(21), -1./7*sqrt(21), -1, 1}; 
			_w = {32./45, 49./90, 49./90, .1, .1}; 
			break; 				
		}

		case 6: {
			x = {-1., -sqrt(1./21*(7.+2*sqrt(7))), -sqrt(1./21*(7-2*sqrt(7))), 
				sqrt(1./21*(7-2*sqrt(7))), sqrt(1./21*(7.+2*sqrt(7))), 1.}; 
			_w = {1./15, 1./30*(14-sqrt(7)), 1./30*(14 + sqrt(7)), 
					1./30*(14 + sqrt(7)), 1./30*(14-sqrt(7)), 1./15}; 
			break; 
		}

		default : {
			ERROR("lobatto integration of order " << _p << " not implemented"); 
		}
	}

	for (int i=0; i<x.GetSize(); i++) {
		_x.Append(Point({x[i]})); 
		// _w[i] *= .5; 
	}
}

QuadTensorProduct::QuadTensorProduct(Quadrature *q, int dim) 
	: Quadrature(q->GetOrder()) {
	_dim = dim; 
	_q = q; 

	if (_dim == 1) {
		for (int i=0; i<_q->NumPoints(); i++) {
			_x.Append(_q->X(i)); 
			_w.Append(_q->Weight(i)); 
		}
	} else if (_dim == 2) {
		for (int i=0; i<_q->NumPoints(); i++) {
			for (int j=0; j<_q->NumPoints(); j++) {
				_x.Append({_q->X(i)[0], _q->X(j)[0]}); 
				_w.Append(_q->Weight(i) * _q->Weight(j)); 
			}
		}
	} else if (_dim == 3) {
		for (int i=0; i<_q->NumPoints(); i++) {
			for (int j=0; j<_q->NumPoints(); j++) {
				for (int k=0; k<_q->NumPoints(); k++) {
					_x.Append({_q->X(i)[0], _q->X(j)[0], _q->X(k)[0]}); 
					_w.Append(_q->Weight(i) * _q->Weight(j) * _q->Weight(k)); 
				}
			}
		}
	} else {
		ERROR("QuadTensorProduct not implemented for dim = " << _dim); 
	}
}

QuadTri::QuadTri(int p) : Quadrature(p) {
	switch(p) {
		case 1: {
			_x.Append({1./3, 1./3}); 
			_w.Append(.5); 
			break; 
		}

		case 4: {
			_x = {
				{0.10810301816807022736,  0.44594849091596488632}, 
				{0.44594849091596488632,  0.10810301816807022736}, 
				{0.44594849091596488632,  0.44594849091596488632}, 
				{0.81684757298045851308,  0.091576213509770743460}, 
				{0.091576213509770743460,  0.81684757298045851308}, 
				{0.091576213509770743460,  0.091576213509770743460}
			}; 

			_w = {
				0.22338158967801146570, 
				0.22338158967801146570, 
				0.22338158967801146570, 
				0.10995174365532186764, 
				0.10995174365532186764, 
				0.10995174365532186764
			}; 
			for (int i=0; i<_w.GetSize(); i++) {
				_w[i] /= 2.; 
			}
			break; 
		}

		case 6: {
			_x = {
				{0.87382197101699554332,   0.063089014491502228340}, 
				{0.063089014491502228340,  0.873821971016995543320}, 
				{0.063089014491502228340,  0.063089014491502228340}, 
				{0.50142650965817915742,   0.249286745170910421290}, 
				{0.24928674517091042129,   0.501426509658179157420}, 
				{0.24928674517091042129,   0.249286745170910421290}, 
				{0.053145049844816947353,  0.310352451033784405420}, 
				{0.31035245103378440542,   0.053145049844816947353}, 
				{0.053145049844816947353,  0.636502499121398647230}, 
				{0.31035245103378440542,   0.636502499121398647230}, 
				{0.63650249912139864723,   0.053145049844816947353}, 
				{0.63650249912139864723,   0.310352451033784405420}
			}; 

			_w = {
				0.050844906370206816921,
				0.050844906370206816921,
				0.050844906370206816921,
				0.11678627572637936603,
				0.11678627572637936603,
				0.11678627572637936603,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194
			}; 
			for (int i=0; i<_w.GetSize(); i++) {
				_w[i] /= 2.; 
			}
			break; 
		}

		default: {
			WARNING("triangular quadrature order " << p << " not implemented. defaulting to order 6."); 
			_x = {
				{0.87382197101699554332,   0.063089014491502228340}, 
				{0.063089014491502228340,  0.873821971016995543320}, 
				{0.063089014491502228340,  0.063089014491502228340}, 
				{0.50142650965817915742,   0.249286745170910421290}, 
				{0.24928674517091042129,   0.501426509658179157420}, 
				{0.24928674517091042129,   0.249286745170910421290}, 
				{0.053145049844816947353,  0.310352451033784405420}, 
				{0.31035245103378440542,   0.053145049844816947353}, 
				{0.053145049844816947353,  0.636502499121398647230}, 
				{0.31035245103378440542,   0.636502499121398647230}, 
				{0.63650249912139864723,   0.053145049844816947353}, 
				{0.63650249912139864723,   0.310352451033784405420}
			}; 

			_w = {
				0.050844906370206816921,
				0.050844906370206816921,
				0.050844906370206816921,
				0.11678627572637936603,
				0.11678627572637936603,
				0.11678627572637936603,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194,
				0.082851075618373575194
			};
			for (int i=0; i<_w.GetSize(); i++) {
				_w[i] /= 2.; 
			}
		}
	}

	CHECK(_x.GetSize() == _w.GetSize()); 
}

QuadRules::QuadRules() {
	_leg.Resize(DIM+1);
	_lob.Resize(DIM+1); 
	for (int i=0; i<_leg.GetSize(); i++) {
		_leg[i].Resize(QUAD_ARRAY_SIZE); 
		_leg[i] = NULL; 

		_lob[i].Resize(QUAD_ARRAY_SIZE); 
		_lob[i] = NULL; 
	}
	_tri.Resize(QUAD_ARRAY_SIZE); 
	_tri = NULL; 
}

QuadRules::~QuadRules() {
	for (int i=0; i<DIM; i++) {
		for (int j=0; j<_leg[i].GetSize(); j++) {
			if (_leg[i][j]) delete _leg[i][j]; 
		}
	}

	for (int i=0; i<DIM; i++) {
		for (int j=0; j<_lob[i].GetSize(); j++) {
			if (_lob[i][j]) delete _lob[i][j]; 
		}
	}

	for (int i=0; i<_tri.GetSize(); i++) {
		if (_tri[i]) delete _tri[i]; 
	}
}

Quadrature* QuadRules::Get(int geom, int order, int type) {
	if (order > QUAD_ARRAY_SIZE) {
		ERROR("quadrature order " << order << " exceeds quadrature array size of " << QUAD_ARRAY_SIZE); 
	}
	if (geom == LINE || geom == QUAD || geom == HEX) {
		int dim; 
		if (geom == LINE) dim = 1; 
		else if (geom == QUAD) dim = 2; 
		else if (geom == HEX) dim = 3; 

		if (type == LEGENDRE) {
			if (!_leg[dim][order]) {
				_leg[dim][order] = new QuadTensorProduct(GetLegendre(order), dim); 
			}
			return _leg[dim][order]; 
		} else if (type == LOBATTO) {
			if (!_lob[dim][order]) {
				_lob[dim][order] = new QuadTensorProduct(GetLobatto(order), dim); 
			}
			return _lob[dim][order]; 
		} else {
			ERROR("type " << type << " not defined"); 
		}
	} else if (geom == TRI) {
		if (!_tri[order]) {
			_tri[order] = new QuadTri(order); 
		}
		return _tri[order]; 
	} else {
		ERROR("geom " << geom << " not defined"); 
	}

	// avoids warning 
	Quadrature* Quadrature; 
	return Quadrature; 
}

Quadrature* QuadRules::GetLegendre(int order) {
	if (!_leg[0][order]) {
		_leg[0][order] = new QuadLegendre(order); 
	}
	return _leg[0][order]; 
}

Quadrature* QuadRules::GetLobatto(int order) {
	if (!_lob[0][order]) {
		_lob[0][order] = new QuadLobatto(order); 
	}
	return _lob[0][order]; 
}

} // end namespace fem 