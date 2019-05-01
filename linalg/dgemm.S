.text
.align 2

#include "rvv.h"

.globl MatMult_RV
.type  MatMult_RV,@function

# N(a0), M(a1), K(a2), A(a3), B(a4), C(a5) 

MatMult_RV: 
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		SCALAR | FP | W64
	)
	slli t1, a1, 3 # stride in A 
	slli t2, a2, 3 # stride in B 
	mul a6, t2, a1 # 8*K * M 
	addi a6, a6, -8 # 8*K*(M-1) 
	add a7, a4, x0 # copy B 
batch: 
	add t5, a2, x0 # copy K 	
rows:
	setvl(t0, a0) # set vl to number of rows 
	vslide v0, v0, t0 # zero v0 
	add t4, a1, x0 # copy M 
loop: 
	vld v3, 0(a4) # load entry of B 
	vlds v1, 0(a3), t1 # column of A 
	vmadd v0, v1, v3, v0 # accumulate A*B 
	add a4, a4, t2 # increment B 
	addi a3, a3, 8 # increment A 
	addi t4, t4, -1 # decrement cols of A 
	bnez t4, loop 
	vsts v0, 0(a5), t2 # store to C 
	sub a3, a3, t1 # reset A
	sub a4, a4, a6 # reset B to next column 
	addi a5, a5, 8 # increment C 
	addi t5, t5, -1 # decrement K 
	bnez t5, rows 
	sub a0, a0, t0 # decrement N by vl 
	add a4, a7, x0 # reset B to start 
	mul t6, t0, t1 # (vl-1)*8M 
	add a3, a3, t6 # update A 
	addi t0, t0, -1 # vl - 1
	mul t6, t0, t2 # (vl-1)*8K 
	add a5, a5, t6 # update C 
	bnez a0, batch
ret
