#include "Array.hpp"

namespace fem 
{

template<> 
void Array<double>::Resize(int N) {
	_vector.resize(N);
	for (int i=0; i<N; i++) {
		_vector[i] = 0.; 
	} 
}

template<>
void Array<int>::Resize(int N) {
	_vector.resize(N,0); 
	for (int i=0; i<N; i++) {
		_vector[i] = 0; 
	} 
}

} // end namespace fem 
