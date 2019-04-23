.text
.align 2

#include "rvv.h"

.globl AddFromDofs_RV
.type  AddFromDofs_RV,@function

# N(0), dofs(a1), v(a2), out(a3) 

AddFromDofs_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | INT | W32
	)
loop:
	setvl(t0, a0) 
	vld v0, 0(a2) 
	vld v3, 0(a1) 
	vsli v3, v3, 3
	vldx v1, 0(a3), v3 
	vadd v0, v0, v1 
	vstx v0, 0(a3), v3
	slli t1, t0, 3 
	slli t2, t0, 2
	add a1, a1, t2 
	add a2, a2, t1 
	sub a0, a0, t0 
	bnez a0, loop 
ret 
