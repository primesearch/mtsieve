/* avx_get.S -- (C) Mark Rodenkirch, June 2018

   All AVX functions require parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[16];

   void  avx_get_all(double *regs)
      This copies the 64 doubles contained in the 16 ymm registers.
      
   void  avx_get_16a(double *regs)
      This copies the 16 doubles contained in ymm12-ymm15, which contain
      the results of the call to avx_mulmod.

   Many thanks to Ernst Meyer for this code.
         
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/
  .text
  
#ifdef _WIN64
#define ARG1      %rcx
#else
#define ARG1      %rdi
#endif

  .p2align 4,,15
  
  .globl _avx_get_all
  .globl avx_get_all

_avx_get_all:
avx_get_all:

   vmovapd	      %ymm0 , 0x0000(ARG1)        
   vmovapd	      %ymm1 , 0x0020(ARG1)        
   vmovapd	      %ymm2 , 0x0040(ARG1)        
   vmovapd	      %ymm3 , 0x0060(ARG1)
   
   vmovapd	      %ymm4 , 0x0080(ARG1)        
   vmovapd	      %ymm5 , 0x00a0(ARG1)        
   vmovapd	      %ymm6 , 0x00c0(ARG1)        
   vmovapd	      %ymm7 , 0x00e0(ARG1)
   
   vmovapd	      %ymm8 , 0x0100(ARG1)          
   vmovapd	      %ymm9 , 0x0120(ARG1)          
   vmovapd	      %ymm10, 0x0140(ARG1)          
   vmovapd	      %ymm11, 0x0160(ARG1)
   
   vmovapd	      %ymm12, 0x0180(ARG1)          
   vmovapd	      %ymm13, 0x01a0(ARG1)          
   vmovapd	      %ymm14, 0x01c0(ARG1)          
   vmovapd	      %ymm15, 0x01e0(ARG1)
   
   ret
   
  .globl _avx_get_16a
  .globl avx_get_16a

_avx_get_16a:
avx_get_16a:
   
   vmovapd	      %ymm12, 0x0000(ARG1)          
   vmovapd	      %ymm13, 0x0020(ARG1)          
   vmovapd	      %ymm14, 0x0040(ARG1)          
   vmovapd	      %ymm15, 0x0060(ARG1)
   
   ret
   