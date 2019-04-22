#ifndef __RVV_H
#define __RVV_H

// vtypes
#define INT     (0 << 6) // integer
#define UINT    (1 << 6) // unsigned integer
#define FP      (3 << 6) // floating point
#define SCALAR  (0 << 11)
#define VECTOR  (4 << 11)
#define W128    48 // 128 bits 
#define W64     32 // 64 bits
#define W32     24 // 32 bits
#define W16     16 // 16 bits 
#define W8      8  // 8 bits

#define setvcfg(vcfg, vtype0, vtype1, vtype2, vtype3) \
  li t0, ((vtype0) | ((vtype1) << 16) | ((vtype2) << 32) | ((vtype3) << 48)) ; \
  csrw vcfg, t0 

#define setvl(rd, rs) \
  csrw vl, rs; \
  csrr rd, vl

#endif // __RVV_H
