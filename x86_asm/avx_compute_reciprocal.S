/* avx_compute_reciprocal.S -- (C) Mark Rodenkirch, June 2018

   All AVX functions require parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[16];
      
   void  avx_compute_reciprocal(double *p, double *reciprocal)
      Computes 1/p for 16 distinct p and stores in the reciprocal array.
      The reciprocal is needed for calls to avx_mulmod and avx_powmod.

   Many thanks to Ernst Meyer for this code.
         
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
  
  .globl _avx_compute_reciprocal
  .globl avx_compute_reciprocal

_avx_compute_reciprocal:
avx_compute_reciprocal:

   mov            $2, %r10
   mov            $3, %r11
   
   cvtsi2sd       %r10, %xmm5
   cvtsi2sd       %r11, %xmm6

   // move 2.0 and 3.0 into ymm14 and ymm15
   vbroadcastsd	%xmm5, %ymm14
   vbroadcastsd	%xmm6, %ymm15
   
   // load inputs-to-be-inverted into 4 AVX registers
   vmovapd	      0x00(ARG1), %ymm4   
   vmovapd	      0x20(ARG1), %ymm5 
   vmovapd	      0x40(ARG1), %ymm6   
   vmovapd	      0x60(ARG1), %ymm7 
   
   // convert to SP
   vcvtpd2ps	   %ymm4, %xmm0	    
   vcvtpd2ps	   %ymm5, %xmm1     
   vcvtpd2ps	   %ymm6, %xmm2     
   vcvtpd2ps	   %ymm7, %xmm3     
   
   // ainv := approx 1/p to 11-12 bits of precision
   vrcpps		   %xmm0, %xmm0	    
   vrcpps		   %xmm1, %xmm1     
   vrcpps		   %xmm2, %xmm2     
   vrcpps		   %xmm3, %xmm3     
   
   // convert ~1/p back to DP
   vcvtps2pd	   %xmm0, %ymm0	    
   vcvtps2pd	   %xmm1, %ymm1     
   vcvtps2pd	   %xmm2, %ymm2     
   vcvtps2pd	   %xmm3, %ymm3     
   
   // 1st NR iteration gives ~23 bits of precision
   
   // make a copy of ainv
   vmovapd	      %ymm0, %ymm8 	    
   vmovapd	      %ymm1, %ymm9        
   vmovapd	      %ymm2, %ymm10       
   vmovapd	      %ymm3, %ymm11       
   
   // 2 - p*ainv, overwrites ainv
   vfnmadd132pd	%ymm4, %ymm14, %ymm0
   vfnmadd132pd	%ymm5, %ymm14, %ymm1
   vfnmadd132pd	%ymm6, %ymm14, %ymm2
   vfnmadd132pd	%ymm7, %ymm14, %ymm3

   // ainv*(2 - p*ainv) = 1/p accurate to ~23 bits
   vmulpd		   %ymm0, %ymm8, %ymm0   
   vmulpd		   %ymm1, %ymm9, %ymm1   
   vmulpd		   %ymm2, %ymm10, %ymm2   
   vmulpd		   %ymm3, %ymm11, %ymm3   
   
   // 3rd-order update of 23-bit result needs just 2 FMA, 1 SUB, 1 MUL:
   
   // make a copy of ainv
   vmovapd	      %ymm0, %ymm8        
   vmovapd	      %ymm1, %ymm9        
   vmovapd	      %ymm2, %ymm10       
   vmovapd	      %ymm3, %ymm11       
   
   // 1st FMA overwrites d data (inputs) with (3 - d*ainv)
   vfnmadd132pd	%ymm0, %ymm15, %ymm4
   vfnmadd132pd	%ymm1, %ymm15, %ymm5
   vfnmadd132pd	%ymm2, %ymm15, %ymm6
   vfnmadd132pd	%ymm3, %ymm15, %ymm7
   
   // Subtract 3 from (3 - p*ainv) to get -y = -p*ainv terms in y0-3
   vsubpd		   %ymm15, %ymm4, %ymm0   
   vsubpd		   %ymm15, %ymm5, %ymm1   
   vsubpd		   %ymm15, %ymm6, %ymm2   
   vsubpd		   %ymm15, %ymm7, %ymm3   
   
   // Positive-product FMA gives (3 - y*(3 - p*ainv)) in y0-3
   vfmadd132pd	   %ymm4, %ymm15, %ymm0
   vfmadd132pd	   %ymm5, %ymm15, %ymm1
   vfmadd132pd	   %ymm6, %ymm15, %ymm2
   vfmadd132pd	   %ymm7, %ymm15, %ymm3
   
   // ainv*(3 - y*(3 - d*ainv)) = 1/p accurate to ~53 bits
   vmulpd		   %ymm0, %ymm8, %ymm8   
   vmulpd		   %ymm1, %ymm9, %ymm9   
   vmulpd		   %ymm2, %ymm10, %ymm10  
   vmulpd		   %ymm3, %ymm11, %ymm11  
   
   // Write the reciprocals to memory, will need these not just now but once more in the final error-correction step:
   vmovapd        %ymm8,  0x00(ARG2)
	vmovapd        %ymm9,  0x20(ARG2)
	vmovapd        %ymm10, 0x40(ARG2)
	vmovapd        %ymm11, 0x60(ARG2)
   
   ret
   