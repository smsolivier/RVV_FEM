.text
.align 2

.globl VectorAdd_RV
.type  VectorAdd_RV,@function

# N(0), a(1), b(2) 

VectorAdd_RV:
	slli a0, a0, 3 # length of vector in bytes 
	add t0, a0, a1 # end pointer for a 
	# add vl, x0, a0 
	# vld v0, 0(a1) 
	# vld v1, 0(a2) 
	# vadd v2, v0, v1 
	# vstx v3, 0, a1
loop:
	fld ft0, 0(a1) 
	fld ft1, 0(a2) 
	fadd.d ft0, ft0, ft1 
	fsd ft0, 0(a1) 
	addi a1, a1, 8 
	addi a2, a2, 8
	bne a1, t0, loop 
ret 
