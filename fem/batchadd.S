.text
.align 2

#include "rvv.h"

.globl BatchAdd_RV
.type  BatchAdd_RV,@function

# N(a0), B(a1), dofs(a2), ball(a3), b(a4)  

BatchAdd_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | INT | W32
	)
	addi t1, a0, 0 # save N 
	slli t2, a0, 2 # 4N 
	slli t3, a0, 3 # 8N 
batch: 
	addi a0, t1, 0 # reset N 
loop:
	setvl(t0, a0) # set vl to dofs/element 
	vld v3, 0(a2) # load dofs 
	vsli v3, v3, 3
	vld v0, 0(a3) # load ball 
	vldx v1, 0(a4), v3 # load b 
	vadd v2, v0, v1 # add ball and b 
	vstx v2, 0(a4), v3 # store back to b
	sub a0, a0, t0 # decrement N by vl 
	slli t4, t0, 2 # 4vl 
	add a2, a2, t4 # update dofs 
	slli t4, t0, 3 # 8vl 
	add a3, a3, t4 # update ball 
	bnez a0, loop 
	addi a1, a1, -1 # decrement batches 
	bnez a1, batch 
ret 