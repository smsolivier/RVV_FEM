.text
.align 2

#include "rvv.h"

.globl MatVec_RV
.type  MatVec_RV,@function

# N(0), M(1), A(2), x(3), b(4)

MatVec_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		SCALAR | FP | W64)
	slli t0, a1, 3 # offset for A 
rows: 
	setvl(t1, a0)
	vinsert v3, x0, x0 
	vadd v0, v3, v3 # zero vector 
	add t2, a2, x0 
	add t3, a3, x0 
	add t4, a1, x0 
loop:
	vlds v1, 0(t2), t0 # col of A 
	vld v3, 0(t3) # entry of x 
	vmadd v0, v1, v3, v0 
	addi t2, t2, 8 
	addi t3, t3, 8
	addi t4, t4, -1 
	bnez t4, loop 
	vst v0, 0(a4)
	mul t3, t1, t0 
	add a2, a2, t3 
	slli t3, t1, 3 
	add a4, a4, t3 
	sub a0, a0, t1 
	bnez a0, rows 
ret
