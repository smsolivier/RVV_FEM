#include "Point.hpp"

namespace fem 
{

Point::Point(std::initializer_list<double> list) {
	int ind = 0; 
	for (auto i=list.begin(); i != list.end(); i++) {
		CHECK(ind < DIM); 
		_x[ind] = *i; 
		ind++; 
	}

	// fill in rest of values to be zero 
	while (ind < DIM) {
		_x[ind] = 0.; 
		ind++; 
	}
}

std::ostream& operator<<(std::ostream& out, const Point& p) {
	int w = 9; 
	out << "(" << std::setw(w) << p[0] << ", " << std::setw(w) <<
		p[1] << ", " << std::setw(w) << p[2] << ")"; 
	return out; 
}

} // end namespace fem 