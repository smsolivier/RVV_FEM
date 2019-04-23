.text
.align 2

#include "rvv.h"

.globl VectorSub2_RV
.type  VectorSub2_RV,@function

# N(0), a(1), b(2) 

VectorSub2_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | FP | W64
	)
loop:
	setvl(t0, a0) 
	vld v0, 0(a1) 
	vld v1, 0(a2) 
	vsub v0, v0, v1 
	vst v0, 0(a3) 
	slli t1, t0, 3 
	add a1, a1, t1 
	add a2, a2, t1 
	add a3, a3, t1 
	sub a0, a0, t0 
	bnez a0, loop 
ret 
