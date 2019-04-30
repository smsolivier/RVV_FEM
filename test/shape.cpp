#include "FEM.hpp"

using namespace std; 
using namespace fem; 

Point x = {0.25, 0.1}; 
bool Shape(Element& el) {
	Array<double> ans = {0.00421875, -0.00703125, 0.00859375, -0.00515625, 
		-0.0421875, 0.154688, 0.0515625, -0.0928125, 0.928125}; 

	Vector shape; 
	el.CalcShape(x, shape); 

	for (int i=0; i<shape.GetSize(); i++) {
		if (fabs(shape[i] - ans[i]) > 1e-6) {
			return false; 
		}
	}
	return true; 
}

bool GradShape(Element& el) {
	Array<double> ans = {0.01125, -0.03375, 0.04125, -0.01375, 0.0225, 0.7425, -0.0275, -0.2475, -0.495, 
0.0375, -0.0625, 0.09375, -0.05625, -0.375, -0.03125, 0.5625, 0.01875, -0.1875}; 
	
	Matrix gshape; 
	el.CalcGradShape(x, gshape); 

	for (int i=0; i<gshape.GetSize(); i++) {
		if (fabs(gshape[i] - ans[i]) > 1e-6) {
			return false; 
		}
	}
	return true; 
}

int main() {
	MeshNode n0 = {0, {0,0}, fem::INTERIOR}; 
	MeshNode n1 = {1, {2,0}, fem::INTERIOR}; 
	MeshNode n2 = {2, {2,1}, fem::INTERIOR}; 
	MeshNode n3 = {3, {0,1}, fem::INTERIOR};
	LagrangeQuad el({n0, n1, n2, n3}, 2);

	TEST(Shape(el), "calc shape"); 
	TEST(GradShape(el), "grad shape"); 
}