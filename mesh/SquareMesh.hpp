#pragma once 

#include "General.hpp"
#include "Mesh.hpp"
#ifdef USE_MPI 
#include <mpi.h>
#endif

namespace fem 
{

/// create a square mesh 
class SquareMesh : public Mesh {
public:
	/// default constructor 
	SquareMesh() { } 
	/// constructor 
	/** \param Nx number of elements in \f$ x \f$ 
		\param Ny number of elements in \f$ y \f$ 
		\param corner upper right corner location (assumes lower left is origin) 
	*/ 
	SquareMesh(int Nx, int Ny, Point corner={1,1}) {SquareMesh(Nx, Ny, {0,0}, corner); }
	/// constructor but specify low corner as well 
	SquareMesh(int Nx, int Ny, Point low, Point high, 
		Array<int> bcs={DIRICHLET, DIRICHLET, DIRICHLET, DIRICHLET}); 

	/// returns number of elements in x and y 
	int GetNx() const {return _Nx; }
	int GetNy() const {return _Ny; }
	/// returns the low corner 
	Point GetLowCorner() const {return _low; }
	/// returns the high corner 
	Point GetHighCorner() const {return _high; }
private: 
	/// number of elements in x 
	int _Nx; 
	/// number of elements in y 
	int _Ny; 
	/// coordinates of low corner 
	Point _low; 
	/// coordinates of high corner 
	Point _high; 
}; 

} // end namespace fem 