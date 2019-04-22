.text
.align 2

.globl MatVec_RV
.type  MatVec_RV,@function

# N(0), M(1), A(2), x(3), b(4)

# MatVec_RV:
# 	slli t0, a0, 3 # height of A in bytes 
# 	slli t1, a1, 3 # width of A in bytes 
# 	add t2, a3, t1 
# 	add t3, a4, t0 
# row: 
# 	fld ft0, 0(a4) 
# dot:
# 	fld ft1, 0(a2) 
# 	fld ft2, 0(a3) 
# 	fmadd.d ft0, ft1, ft2, ft0 
# 	addi a2, a2, 8 
# 	addi a3, a3, 8 
# 	bne a3, t2, dot 
# 	fsd ft0, 0(a4) 
# 	addi a4, a4, 8
# 	sub a3, a3, t1
# 	bne a4, t3, row 
# ret 
MatVec_RV:
	