#include "FEM.hpp"

using namespace std; 
using namespace fem; 

Point x = {.25, .25}; 

int main(int argc, char* argv[]) {
	int p = 1; 
	if (argc > 1) p = atoi(argv[1]); 
	MeshNode n0 = {0, {0,0}, fem::INTERIOR}; 
	MeshNode n1 = {1, {2,0}, fem::INTERIOR}; 
	MeshNode n2 = {2, {2,1}, fem::INTERIOR}; 
	MeshNode n3 = {3, {0,1}, fem::INTERIOR};
	LagrangeQuad el({n0, n1, n2, n3}, p);

	Vector shape; 
	Matrix gshape; 
	HWCounter hwc_s; 
	for (int i=0; i<100; i++) {
		el.CalcShape(x, shape); 		
	}
	hwc_s.Read(); 

	HWCounter hwc_g; 
	for (int i=0; i<100; i++) {
		el.CalcGradShape(x, gshape); 		
	}
	hwc_g.Read(); 

	hwc_s.PrintStats("calc shape"); 
	hwc_g.PrintStats("calc grad shape"); 
}