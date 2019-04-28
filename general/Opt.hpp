#pragma once 

#define VECTORIZATION
// #define RV_SMATVEC

#ifdef VECTORIZATION 
// matrix optimizations 
#define RV_MATVEC 
// #define RV_MATMULT // not implemented 

// vector optimizations 
#define RV_VECADD 
#define RV_VECSUB 
#define RV_VECDIV 
#define RV_VECMUL 
#define RV_VECSCALE 
#define RV_VECOP 
#define RV_VECDOT 
#define RV_GETDOFS 
#define RV_ADDDOFS

// fematrix optimizations 
#define RV_MVOUTER

#endif