#pragma once 

#include "General.hpp"
#include "GridFunction.hpp"
#include "Array.hpp"
#ifdef USE_MPI 
#include <mpi.h>
#endif

namespace fem 
{

/// write to vtk 
class Writer {
public:
	/// constructor. provide base name for output files 
	Writer(std::string name="solution"); 
	/// add a solution variable to the output list
	/** \param gf GridFunction memory location 
		\param name name to write out 
	*/  
	void Add(GridFunction& gf, std::string name); 
	/// set the write frequency 
	/** write every f calls to Write */ 
	void SetFreq(int f) {_f = f; }

	/// write all variables to gmsh 
	/** \param force force output regardless of write frequency settings */ 
	void WriteGMSH(bool force=false); 

	/// write all variables to VTK 
	void Write(bool force=false); 
protected:
	/// output frequency. write every _f times 
	int _f; 
	/// store pointers to GridFunction's 
	Array<GridFunction*> _gf; 
	/// store the corresponding output names 
	Array<std::string> _names; 
	/// store the base name for outputting 
	std::string _base_name; 
	/// number of times Write has been called 
	int _count; 
	/// number of files written 
	int _writes; 
}; 

} // end namespace fem 