// #pragma once 
#include <iostream> 
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector> 
#include <array> 
#include <string> 
#include <chrono> 
#include <cstring>
#ifdef USE_MPI 
#include <mpi.h>
#endif

#define DIM 3 
#define MESH_DIR ""
#define INTEGRATION_ORDER 8
#define INTEGRATION_TYPE 0

// exit with error and backtracing 
#define EXIT {\
	std::cout << "encountered FATAL ERROR... aborting" << std::endl; \
	exit(1); \
}

#define EXIT_NOBACKTRACE {\
	std::cout << "encountered FATAL ERROR... aborting" << std::endl; \
	exit(1); \
}

// give an error message and exit 
#define ERROR(msg) {\
	std::cout << "FEM ERROR: \n --> " << msg << "\n ... in " << __PRETTY_FUNCTION__ \
	<< "\n ... from " << __FILE__ << " " << __LINE__ << std::endl; \
	EXIT; \
}

#define ERRORCHECK(chk, msg) {\
	if ((bool)(chk)) {\
		std::cout << "ERROR: " << msg << " (" << #chk << ")" << \
		"\n\tin " << __PRETTY_FUNCTION__ \
		<< "\n\tfrom " << __FILE__ << " " << __LINE__ << std::endl; \
		EXIT; \
	}\
}

// output a warning 
#ifdef USE_WARNINGS 
#define WARNING(msg) {\
	std::cout << "WARNING:\n --> " << msg << "\n ... in " << __PRETTY_FUNCTION__ \
	<< "\n ... from " << __FILE__ << " " << __LINE__ << std:: endl; \
}
#else 
#define WARNING(msg) 
#endif

#define ASSERT(check) {\
	if (!(bool)(check)) {\
		std::cout << "ASSERTION FAILED: " << #check << " is false" << std::endl \
		<< "\tin " << __PRETTY_FUNCTION__ << "\n\tfrom " << __FILE__ \
		<< " " << __LINE__ << std::endl; \
		EXIT; \
	}\
}

#ifndef NDEBUG 
#define CHECK(chk) {\
	if (!(bool)(chk)) {\
		std::cout << "CHECK FAILED: " << #chk << " is false" << std::endl \
		<< "\tin " << __PRETTY_FUNCTION__ << "\n\tfrom " << __FILE__ \
		<< " " << __LINE__ << std::endl; \
		EXIT; \
	}\
}
#else 
#define CHECK(chk) 
#endif

#ifndef NDEBUG 
#define CHECKMSG(chk, msg) {\
	if (!(bool)(chk)) {\
		std::cout << "CHECK FAILED: " << #chk << " is false" \
		<< std::endl << " --> " << msg << std:: endl \
		<< " ... in " << __PRETTY_FUNCTION__ << "\n ... from " << __FILE__ \
		<< " " << __LINE__ << std::endl; \
		EXIT; \
	}\
}
#else 
#define CHECKMSG(chk, msg)
#endif

#ifdef TESTING
#define TEST(chk, msg) {\
	std::cout << msg << ": "; \
	if ((bool)(chk)) {\
		std::cout << "pass" << std::endl; \
	} else {\
		std::cout << "fail" << std:: endl; \
		EXIT_NOBACKTRACE; \
	}\
}
#else 
#define TEST(chk, msg)
#endif

#define EQUAL(a, b) (fabs(a-b) < 1e-7)

#ifndef NDEBUG 
#define PRINTVAR(a) {std::cout << #a << " = " << a << std::endl; }
#else 
#define PRINTVAR(a) 
#endif