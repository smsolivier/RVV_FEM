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
batch: 
	setvl(t0, a1) # set according to batch size 
	mul t5, t2, t0 # amount to skip in dofs array 
	mul a6, t1, t0 # amount to skip in mats array 
	vlds v3, 0(a3), t2 # load dofs
	vsli v4, v3, 3 # dofs for x (must be reset) 
	vsli v5, v3, 3 # dofs for b (incremented with rows)
	vsli v3, v3, 3 # base copy 
	add t6, a2, x0 # save location of mats 
rows:
	add t3, t4, x0 # copy of N 
	vslide v2, v2, t0 # zero v2 
loop:
	vlds v0, 0(t6), t1 # first element of matrices 
	vldx v1, 0(a4), v4 # load first element of x
	vmadd v2, v0, v1, v2 # multiply Ax
	vaddi v4, v4, 8 # increment x locations 
	addi t3, t3, -1 # decrement N 
	addi t6, t6, 8 # move to next column of A 
	bnez t3, loop
	vldx v0, 0(a5), v5
	vadd v2, v0, v2 
	vstx v2, 0(a5), v5 # store to b 
	addi a0, a0, -1 # decrement rows 
	vaddi v4, v3, 0 # reset x dofs 
	vaddi v5, v5, 8
	bnez a0, rows
	sub a1, a1, t0 # decrement by vl  
	add a0, t4, x0 # reset N for rows 
	add a3, a3, t5 # increment dofs 
	add a2, a2, a6 # increase mats array by vl 
	bnez a1, batch
ret	