.text
.align 2

#include "rvv.h"

.globl VectorScale_RV
.type  VectorScale_RV,@function

# N(0), v(1), alpha(fa0) 

# VectorScale_RV:
# 	slli t0, a0, 3
# loop:
# 	fld f0, 0(a1) 
# 	fmul.d f0, f0, fa0
# 	fsd f0, 0(a1) 
# 	addi a1, a1, 8
# 	addi a0, a0, -1 
# 	bnez a0, loop 
# ret
VectorScale_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		SCALAR | FP | W64
	)
	vld v3, 0(a2) 
loop:
	setvl(t0, a0) 
	vld v0, 0(a1) 
	vmul v0, v0, v3
	vst v0, 0(a1) 
	slli t1, t0, 3 
	add a1, a1, t1 
	sub a0, a0, t0
	bnez a0, loop
ret 
