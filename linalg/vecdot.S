.text
.align 2

#include "rvv.h"

.globl VectorDot_RV
.type  VectorDot_RV,@function

# N(a0), a(a1), b(a2), c(a3) 

VectorDot_RV:
    # configufations                                                                                                  
    # v0-v2: int
    setvcfg(vcfg0,
        VECTOR | FP | W64,                                                                      
        VECTOR | FP | W64,                                                                         
        VECTOR | FP | W64,                                                                         
        SCALAR | FP | W64)
    # Setup: set t0 = MVL by asking for a large value
    li a4, 1024                                                                                                       
    setvl(t0, a4)
    vinsert v3, x0, x0
    # vadd v2, v3, v3 
    setvcfg(vcfg0,
        VECTOR | FP | W64,                                                                      
        VECTOR | FP | W64,                                                                         
        VECTOR | FP | W64,                                                                         
        VECTOR | FP | W64)
#     # Stripmine: accumulate partial results in each element
stripmine:
    setvl(a4, a0)
    # load from A
    vld v0, 0(a1)
    # load from B
    vld v1, 0(a2)
    # do the fused multiply add
    vmadd v2, v0, v1, v2
    # bump pointers and decrement n
    slli a5, a4, 3
    add a1, a1, a5
    add a2, a2, a5
    sub  a0, a0, a4
    # loop until all products are finished
    bnez a0, stripmine

    # Reduction: Do a tree reduction on the full vector length
    # Note: this is clearly not optimal for short vectors, but it works in all cases
    li a5, 1 # Stop condition -> VL = 1
    setvl(a4, t0) # Set VL = MVL
reduction:
    srli a4, a4, 1     # Divide VL by 2
    vslide v1, v2, a4  # Slide the top half of the vector into v1
#     setvl(a4, a4)
#     vadd v2, v1, v2    # Do the add  
#     bne a4, a5, reduction

    # vld v1, 0(a3)      # VL == 1; Get the value stored at C
    # vadd v2, v1, v2    # VL == 1; Add the initial value to the result of the ddot
    # vst v2, 0(a3)      # VL == 1; Store the single element back to memory
ret
