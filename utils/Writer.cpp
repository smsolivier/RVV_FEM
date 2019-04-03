#include "Writer.hpp"
#include "VisitWriter.H"

using namespace std; 

namespace fem 
{

Writer::Writer(string name) {
	_base_name = name; 
	_count = 0; 
	_writes = 0; 
	_f = 1; 
}

void Writer::Add(GridFunction& gf, string name) {
	_gf.Append(&gf); 
	_names.Append(name); 
}

void Writer::WriteGMSH(bool force) {

	if (_count++%_f != 0 && !force) {
		return; 
	}

	const FESpace* space = _gf[0]->GetSpace(); 
	string name = _names[0]; 

	ofstream out(_base_name + to_string(_writes++) + ".msh"); 

	// header 
	out << "$MeshFormat" << endl; 
	out << 2.2 << " " << 0 << " " << 8 << endl; 
	out << "$EndMeshFormat" << endl; 

	// nodes 
	out << "$Nodes" << endl; 
	out << space->GetNumNodes() << endl; 
	for (int i=0; i<space->GetNumNodes(); i++) {
		out << space->GetNode(i).GetGlobalID()+1 << " "; 
		for (int j=0; j<DIM; j++) {
			out << space->GetNode(i).GetX()[j] << " "; 
		}
		if (DIM==2) {
			out << 0; 
		}
		out << endl; 
	}
	out << "$EndNodes" << endl; 

	// elements 
	out << "$Elements" << endl; 
	out << space->GetNumElements() << endl; 
	for (int i=0; i<space->GetNumElements(); i++) {
		const Element& el = space->GetEl(i); 
		out << i+1 << " " << el.GetType()+1 << " " 
			<< 0 << " "; 

		for (int j=0; j<el.GetNumNodes(); j++) {
			out << el[j].GetGlobalID()+1 << " "; 
		}
		out << endl; 
	}
	out << "$EndElements" << endl; 

	// node data 
	out << "$NodeData" << endl; 
	out << 
		"1" << endl << "\"" << name << "\"" << endl << // string tag 
		"1" << endl << "0.0" << endl << // real tag 
		"3" << endl << // number of integer tags 
		_writes-1 << endl << // time step number 
		"1" << endl << // number of components 
		space->GetNumNodes() << endl; // number of node values to write 
	for (int i=0; i<space->GetNumNodes(); i++) {
		out << space->GetNode(i).GetGlobalID()+1 << " " << (*_gf[0])[i] << endl; 
	}
	out << "$EndNodeData" << endl; 

	out.close(); 
}

void Writer::Write(bool force) {
#ifndef USE_RISCV
	if (_count++%_f != 0 && !force) {
		return; 
	}

	const FESpace* space = _gf[0]->GetSpace(); 
	int N = space->GetNumNodes(); 

	int nvars = _gf.GetSize(); 

	// deep copy to float** 
	float** vars = new float*[nvars]; 
	int vardim[nvars]; 
	for (int i=0; i<nvars; i++) {
		const FESpace* fes = _gf[i]->GetSpace(); 
		int dim = fes->GetVDim(); 
		int vdim = (dim > 1) ? 3 : 1; 
		vardim[i] = vdim; 
		vars[i] = new float[vdim*fes->GetNumNodes()]; 
		if (vdim == 1) {
			for (int j=0; j<_gf[i]->GetSize(); j++) {
				vars[i][j] = _gf[i]->operator[](j); 
			} 			
		} else {
			if (dim==2) {
				for (int j=0; j<fes->GetNumNodes(); j++) {
					vars[i][vdim*j] = _gf[i]->operator[](dim*j); 
					vars[i][vdim*j+1] = _gf[i]->operator[](dim*j+1); 
					vars[i][vdim*j+2] = 0; 
				}
			} else {
				for (int j=0; j<fes->GetNumNodes(); j++) {
					for (int d=0; d<3; d++) {
						vars[i][vdim*j+d] = _gf[i]->operator[](dim*j+d); 
					}
				}
			}
		}
	}

	int centering[nvars]; 
	for (int i=0; i<nvars; i++) {
		centering[i] = 1; 
	}

	const char* varnames[nvars]; 
	for (int i=0; i<nvars; i++) {
		varnames[i] = _names[i].c_str(); 
	}

	vector<float> pts(3*N);
	for (int i=0; i<N; i++) {
		Point x = space->GetNode(i).GetX(); 
		int p = 3*space->GetNode(i).GetGlobalID();  
		for (int j=0; j<DIM; j++) {
			pts[p+j] = x[j]; 
		}
	} 

	vector<int> conns; 
	vector<int> cellType; 
	for (int i=0; i<space->GetNumElements(); i++) {
		const Element& el = space->GetEl(i); 

		if (el.GetType() == QUAD) {
			if (space->GetOrder() == 0) {
				ERROR("not supported"); 
			}
			
			else if (space->GetOrder() == 1) {
				for (int j=0; j<4; j++) {
					conns.push_back(el[j].GetGlobalID()); 
				}
				cellType.push_back(VISIT_QUAD); 
			} 

			else if (space->GetOrder() == 2) {
				conns.push_back(el.GetNodeGlobalID(0)); 
				conns.push_back(el.GetNodeGlobalID(4)); 
				conns.push_back(el.GetNodeGlobalID(8)); 
				conns.push_back(el.GetNodeGlobalID(7)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(4)); 
				conns.push_back(el.GetNodeGlobalID(1)); 
				conns.push_back(el.GetNodeGlobalID(5)); 
				conns.push_back(el.GetNodeGlobalID(8)); 
				cellType.push_back(VISIT_QUAD);

				conns.push_back(el.GetNodeGlobalID(8)); 
				conns.push_back(el.GetNodeGlobalID(5)); 
				conns.push_back(el.GetNodeGlobalID(2)); 
				conns.push_back(el.GetNodeGlobalID(6)); 
				cellType.push_back(VISIT_QUAD);

				conns.push_back(el.GetNodeGlobalID(7)); 
				conns.push_back(el.GetNodeGlobalID(8)); 
				conns.push_back(el.GetNodeGlobalID(6)); 
				conns.push_back(el.GetNodeGlobalID(3)); 
				cellType.push_back(VISIT_QUAD); 
			} else if (space->GetOrder() == 3) {
				conns.push_back(el.GetNodeGlobalID(0)); 
				conns.push_back(el.GetNodeGlobalID(4)); 
				conns.push_back(el.GetNodeGlobalID(12)); 
				conns.push_back(el.GetNodeGlobalID(11)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(4)); 
				conns.push_back(el.GetNodeGlobalID(5)); 
				conns.push_back(el.GetNodeGlobalID(13)); 
				conns.push_back(el.GetNodeGlobalID(12)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(5)); 
				conns.push_back(el.GetNodeGlobalID(1)); 
				conns.push_back(el.GetNodeGlobalID(6)); 
				conns.push_back(el.GetNodeGlobalID(13)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(11)); 
				conns.push_back(el.GetNodeGlobalID(12)); 
				conns.push_back(el.GetNodeGlobalID(14)); 
				conns.push_back(el.GetNodeGlobalID(10)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(12)); 
				conns.push_back(el.GetNodeGlobalID(13)); 
				conns.push_back(el.GetNodeGlobalID(15)); 
				conns.push_back(el.GetNodeGlobalID(14)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(13)); 
				conns.push_back(el.GetNodeGlobalID(6)); 
				conns.push_back(el.GetNodeGlobalID(7)); 
				conns.push_back(el.GetNodeGlobalID(15)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(10)); 
				conns.push_back(el.GetNodeGlobalID(14)); 
				conns.push_back(el.GetNodeGlobalID(9)); 
				conns.push_back(el.GetNodeGlobalID(3)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(14)); 
				conns.push_back(el.GetNodeGlobalID(15)); 
				conns.push_back(el.GetNodeGlobalID(8)); 
				conns.push_back(el.GetNodeGlobalID(9)); 
				cellType.push_back(VISIT_QUAD); 

				conns.push_back(el.GetNodeGlobalID(15)); 
				conns.push_back(el.GetNodeGlobalID(7)); 
				conns.push_back(el.GetNodeGlobalID(2)); 
				conns.push_back(el.GetNodeGlobalID(8)); 
				cellType.push_back(VISIT_QUAD); 
			}
		}

		else if (el.GetType() == TRI) {
			for (int i=0; i<3; i++) {
				conns.push_back(el.GetNodeGlobalID(i)); 
			}
			cellType.push_back(VISIT_TRIANGLE); 
		}

		else if (el.GetType() == HEX) {
			for (int i=0; i<8; i++) {
				conns.push_back(el.GetNodeGlobalID(i)); 
			}
			cellType.push_back(VISIT_HEXAHEDRON); 
		}

		else {
			ERROR("element type " << el.GetType() << " not defined"); 
		}
	}

	string fname = _base_name + to_string(_writes++); 
	write_unstructured_mesh(fname.c_str(), 1, N, &(pts[0]), 
		cellType.size(), &(cellType[0]), &(conns[0]), nvars, vardim, centering, 
		varnames, vars); 
#endif
}

} // end namespace fem 