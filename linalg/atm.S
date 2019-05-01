.text
.align 2

#include "rvv.h"

.globl AddTransMult_RV
.type  AddTransMult_RV,@function

# N(a0), M(a1), K(a2), A(a3), B(a4), C(a5) 

AddTransMult_RV: 
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		SCALAR | FP | W64
	)
	slli t1, a1, 3 # stride in A 
	slli t2, a2, 3 # stride in B 
	mul a6, t2, a0 # 8*K*N
	addi a6, a6,-8 # 8*K*(N-1) 
	add a7, a4, x0 # copy of B 
	mul t3, t1, a0 # 8*M*N 
batch: 
	add t5, a2, x0 # copy of K 
rows: 
	setvl(t0, a1) # set to M 
	vslide v0, v0, t0 # zero v0 
	add t4, a0, x0 # copy N 
loop:
	vld v3, 0(a4) # load an entry of B 
	vld v1, 0(a3) # load A 
	vmadd v0, v1, v3, v0 # accumulate A*B 
	add a4, a4, t2 # increment B
	add a3, a3, t1 # next column of A 
	addi t4, t4, -1 # decrement N 
	bnez t4, loop 
	vlds v1, 0(a5), t2 # load C 
	vadd v0, v0, v1 # add to C 
	vsts v0, 0(a5), t2 # store to C 
	sub a4, a4, a6 # reset B to next column 
	sub a3, a3, t3 # reset A to start 
	addi a5, a5, 8 # update C 
	addi t5, t5, -1 # decrement K 
	bnez t5, rows 
	sub a1, a1, t0 # decrement M by vl 
	add a4, a7, x0 # reset B to start 
	slli t6, t0, 3 # 8*vl 
	add a3, a3, t6 # update A 
	addi t0, t0, -1 # vl -1 
	mul t6, t0, t2 # (vl-1)*8K 
	add a5, a5, t6 # update C 
	bnez a1, batch 
ret 