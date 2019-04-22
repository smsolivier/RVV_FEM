.text
.align 2

#include "rvv.h"

.globl VectorOP_RV
.type  VectorOP_RV,@function

# Na(a0), Nb(a1), a(a2), b(a3), M(a4) 

VectorOP_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		SCALAR | FP | W64
	)
aloop:
	add t3, a3, x0 # copy start of b 
	add t2, a1, x0 # copy length of b 
	vld v3, 0(a2) # scalar load a 
bloop: 
	setvl(t0, t2) 
	vld v0, 0(t3) # vector load b 
	vmul v0, v0, v3 
	vst v0, 0(a4) 
	slli t1, t0, 3 
	add a4, a4, t1 
	add t3, t3, t1 
	sub t2, t2, t0 
	bnez t2, bloop 
	addi a0, a0, -1
	addi a2, a2, 8
	bnez a0, aloop 
ret 
