.text
.align 2

#include "rvv.h"

.globl MatMult_RV
.type  MatMult_RV,@function

# N(a0), M(a1), A(a2), B(a3), C(a4) 

MatMult_RV: 
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | FP | W64
	)
ret
