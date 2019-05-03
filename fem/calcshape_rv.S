.text
.align 2

#include "rvv.h"

# Nn(a0), Nc(a1), coef_x(a2), coef_y(a3), shape(a4), x(a5) 

.globl CalcShape_RV
.type  CalcShape_RV,@function

CalcShape_RV: 
	setvcfg(vcfg0, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		VECTOR | FP | W64,
		VECTOR | FP | W64
	)
	setvcfg(vcfg2, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		SCALAR | FP | W64,
		SCALAR | FP | W64
	)
	setvcfg(vcfg4, 
		VECTOR | FP | W64, 
		VECTOR | FP | W64,
		SCALAR | FP | W64,
		SCALAR | FP | W64
	)
	setvl(t0, a0)
	vld v8, 0(a4) # load ones 
	# load x and y in scalar vector registers 
	vld v6, 0(a5) # x 
	vld v7, 8(a5) # y 
	slli t1, a1, 3 # stride in coefs array 
	addi t6, a4, 0 # copy of shape 
batch: 
	setvl(t0, a0) # set to number of shape functions 
	vslide v2, v2, t0 # zero v2 
	vslide v3, v3, t0 # zero v3 
	vslide v4, v8, x0 # reset x to 1 
	vslide v5, v8, x0 # reset y to 1 
	addi t2, a2, 0 # copy of x coefs 
	addi t3, a3, 0 # copy of y coefs
	addi t4, a1, 0 # copy of Nc 
loop:
	vlds v0, 0(t2), t1 # load x coefs 
	vlds v1, 0(t3), t1 # load y coefs 
	vmadd v2, v0, v4, v2 # multiply x*coefs 
	vmadd v3, v1, v5, v3 # multiply y*coefs 
	vmul v4, v4, v6 # accumulate powers of x 
	vmul v5, v5, v7 # accumulate powers of y 
	addi t2, t2, 8 # incremenet x coefs 
	addi t3, t3, 8 # incremenet y coefs 
	addi t4, t4, -1 # decrement number of coefficients 
	bnez t4, loop 
	vmul v2, v2, v3 # multiply shapex and shapey 
	vst v2, 0(a4) 
	sub a0, a0, t0 # decrement Nb by vl 
	mul t5, t0, t1 # vl*Nc*8 
	add a2, a2, t5 # update x coefs 
	add a3, a3, t5 # update y coefs 
	slli t5, t0, 3 # convert vl to bytes 
	add a4, a4, t5 # update shape 
	bnez a0, batch 
ret 
