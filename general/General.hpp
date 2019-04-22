// #pragma once 
#include <iostream> 
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector> 
#include <array> 
#include <string> 
#include <chrono> 
#include <cstring>
#ifdef USE_MPI 
#include <mpi.h>
#endif

#ifndef M_PI 
#define M_PI (3.14159265358979323846) 
#endif

#define DIM 3 
#define MESH_DIR ""
#define INTEGRATION_ORDER 3
#define INTEGRATION_TYPE 0

#define read_csr(reg) ({ unsigned long __tmp; \
  asm volatile ("csrr %0, " #reg : "=r"(__tmp)); \
  __tmp; })

#define write_csr(reg, val) ({ \
  asm volatile ("csrw " #reg ", %0" :: "rK"(val)); })

#define swap_csr(reg, val) ({ unsigned long __tmp; \
  asm volatile ("csrrw %0, " #reg ", %1" : "=r"(__tmp) : "rK"(val)); \
  __tmp; })

#define set_csr(reg, bit) ({ unsigned long __tmp; \
  asm volatile ("csrrs %0, " #reg ", %1" : "=r"(__tmp) : "rK"(bit)); \
  __tmp; })

#define clear_csr(reg, bit) ({ unsigned long __tmp; \
  asm volatile ("csrrc %0, " #reg ", %1" : "=r"(__tmp) : "rK"(bit)); \
  __tmp; })

// exit with error and backtracing 
#define EXIT {\
	std::cout << "encountered FATAL ERROR... aborting" << std::endl; \
	exit(EXIT_FAILURE); \
}

#define EXIT_NOBACKTRACE {\
	std::cout << "encountered FATAL ERROR... aborting" << std::endl; \
	exit(EXIT_FAILURE); \
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