.text
.align 2

#include "rvv.h"

.globl SetEqual_RV
.type  SetEqual_RV,@function

# N(a0), val(a1), data(a2) 

SetEqual_RV:
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		SCALAR | FP | W64
	)
	li t0, 1024
	setvl(t0, t0) 
	vld v3, 0(a1) # load val 
	vslide v0, v0, t0 # zero 
	vadd v0, v0, v3 # set to val 
loop:
	setvl(t0, a0) 
	vst v0, 0(a2) # store to data 
	sub a0, a0, t0 
	slli t0, t0, 3 # convert to double 
	add a2, a2, t0 # update data 
	bnez a0, loop 
ret 
