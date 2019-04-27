.text
.align 2

#include "rvv.h"

.globl MVOuterC_RV
.type  MVOuterC_RV,@function

# N(a0), B(a1), mats(a2), dofs(a3), x(a4), b(a5) 

MVOuterC_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | INT | W32
	)
	setvcfg(vcfg2, 
		VECTOR | INT | W32, 
		VECTOR | INT | W32,
		VECTOR | FP | W64,
		VECTOR | FP | W64
	)
	mul a6, a0, a0 
	slli a6, a6, 3 # 8N*N 
	slli a7, a0, 2 # 4N 
	slli t3, a0, 3 # 8N  
batch: 
	setvl(t0, a1) # set vl to batch size 
	mul t1, a7, t0 # vl * 4N 
	mul t2, a6, t0 # vl * 8N*N 
	addi t5, a0, 0 # rows copy of N
rows:
	vslide v2, v2, t0 # zero v2 
	addi t4, a0, 0 # cols copy of N
loop:
	vlds v3, 0(a3), a7 # load dofs 
	vsli v3, v3, 3 # dofs for x 
	vlds v0, 0(a2), a6 # load element of matrices 
	vldx v1, 0(a4), v3 # load x 
	vmadd v2, v0, v1, v2 # accumulate dot product 
	addi t4, t4, -1 # decrement cols 
	addi a2, a2, 8 # move to next column of mats 
	addi a3, a3, 4 # update dofs 
	bnez t4, loop 
	vsts v2, 0(a5), t3 # store into b 
	addi t5, t5, -1 # decrement rows 
	sub a3, a3, a7 # reset dofs 
	addi a5, a5, 8 # increment b 
	bnez t5, rows 
	sub a1, a1, t0 # decrement batches by vl 
	mul t6, t0, a7 # (vl-1) * 4N 
	add a3, a3, t6 # update dofs 
	addi t0, t0, -1 # subtract one from vl 
	mul t6, t0, a6 # (vl-1) * 8N^2 
	add a2, a2, t6 # update mats to new batch 
	mul t6, t0, t3 # (vl-1) * 8N 
	add a5, a5, t6 # update b 
	bnez a1, batch 
ret 
