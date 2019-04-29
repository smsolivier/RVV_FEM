.text
.align 2

#include "rvv.h"

.globl MVOuter_RV
.type  MVOuter_RV,@function

# N(a0), B(a1), mats(a2), dofs(a3), x(a4), b(a5) 

MVOuter_RV:
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
	mul t1, a0, a0 # N*N 
	slli t1, t1, 3 # N*N*8 
	slli t2, a0, 2 # N*4 
	addi t4, a0, 0 # copy N 
	addi t6, a3, 0 # copy dofs start 
batch: 
	setvl(t0, a1) # set according to batch size 
	addi t6, a3, 0 # copy a3 
rows:
	add a7, t6, x0 # copy dofs location 
	add t3, t4, x0 # copy of N 
	vslide v2, v2, t0 # zero v2 
loop:
	vlds v4, 0(a7), t2 # load next dofs 
	vsli v4, v4, 3 # convert to double 
	vlds v0, 0(a2), t1 # load element of matrices 
	vldx v1, 0(a4), v4 # load element of x
	vmadd v2, v0, v1, v2 # accumulate Ax 
	addi t3, t3, -1 # decrement column number 
	addi a2, a2, 8 # move to next column of A 
	addi a7, a7, 4 # increment dofs 
	bnez t3, loop
	vlds v3, 0(a3), t2 # load b dofs 
	vsli v3, v3, 3 # convert to double
	vldx v0, 0(a5), v3 # load b 
	vadd v2, v0, v2 # add dot product to current b entry 
	vstx v2, 0(a5), v3 # store to b
	addi a3, a3, 4 # update b dofs 
	addi a0, a0, -1 # decrement row number 
	bnez a0, rows
	# do next batch of matvecs 
	sub a1, a1, t0 # decrement by vl  
	add a0, t4, x0 # reset N for rows 
	addi t0, t0, -1 # vl - 1 
	mul t5, t0, t2 # (vl-1)*4N 
	add a3, a3, t5 # increment dofs by (vl-1)*4N 
	mul t5, t0, t1 # (vl-1)*8N^2 
	add a2, a2, t5 # update mats by (vl-1)*8N^2
	bnez a1, batch
ret	