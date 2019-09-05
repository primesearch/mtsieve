/* avx_set_a.S -- (C) Mark Rodenkirch, June 2018

   All AVX functions require parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[16];
      
   void  avx_set_1a(double *a)
      Set all 16 doubles in ymm12-ymm15 to the same a, which must be set
      before calling avx_mulmod or avx_powmod
      
   void  avx_set_16a(double *a)
      Set all 16 doubles in ymm12-ymm15 to 16 distinct a, which must be set
      before calling avx_mulmod or avx_powmod

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

  .text
  
#ifdef _WIN64
#define ARG1      %rcx
#define ARG2      %rdx
#define ARG3      %r8
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#endif

  .p2align 4,,15
  
  .globl _avx_set_1a
  .globl avx_set_1a

_avx_set_1a:
avx_set_1a:
   vbroadcastsd   0(ARG1), %ymm12
   vbroadcastsd   0(ARG1), %ymm13
   vbroadcastsd   0(ARG1), %ymm14
   vbroadcastsd   0(ARG1), %ymm15
   
   ret

  .globl _avx_set_16a
  .globl avx_set_16a

_avx_set_16a:
avx_set_16a:
   vmovapd        0x00(ARG1), %ymm12
   vmovapd        0x20(ARG1), %ymm13
   vmovapd        0x40(ARG1), %ymm14
   vmovapd        0x60(ARG1), %ymm15
   
   
   ret
   