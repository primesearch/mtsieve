/* avx_mulmod.S -- (C) Mark Rodenkirch, May 2018

   All AVX functions require parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[16];

   uint16_t  avx_mulmod(double * p, double *reciprocal)
      Compute a = (a*b)%p for 1 < a,b,p < 2^52 for 16 different a and p, but 1 b.
      
      You must call avx_mulmod_init to compute the reciprocal that is passed
      to this function.

      Upon input ymm4-ymm7 contain the 16 values of b.
      Upon input ymm12-ymm15 contain the 16 values of a.
      Upon output ymm12-ymm15 will contain the 16 remainders of the mulmod.
      Contents of the other ymm registers are lost.
      
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

  .text
  
#ifdef _WIN64
#define ARG1      %rcx
#define ARG2      %rdx
#else
#define ARG1      %rdi
#define ARG2      %rsi
#endif

  .p2align 4,,15
  
  .globl _avx_mulmod
  .globl avx_mulmod

_avx_mulmod:
avx_mulmod:

   // a * b, high 53 bits
   vmulpd         %ymm12, %ymm4, %ymm0
   vmulpd         %ymm13, %ymm5, %ymm1
   vmulpd         %ymm14, %ymm6, %ymm2
   vmulpd         %ymm15, %ymm7, %ymm3
   
   // q = trunc(a * b / p)
   vmulpd         0x00(ARG2), %ymm0, %ymm8    
   vmulpd         0x20(ARG2), %ymm1, %ymm9    
   vmulpd         0x40(ARG2), %ymm2, %ymm10   
   vmulpd         0x60(ARG2), %ymm3, %ymm11   
   
   // a * b, low 53 bits
   vfmsub132pd    %ymm4, %ymm0, %ymm12  
   vfmsub132pd    %ymm5, %ymm1, %ymm13  
   vfmsub132pd    %ymm6, %ymm2, %ymm14  
   vfmsub132pd    %ymm7, %ymm3, %ymm15  
   
   // Round q toward 0 (since these are >= 0 could also round toward -oo)
   vroundpd       $3, %ymm8 , %ymm8       
   vroundpd       $3, %ymm9 , %ymm9       
   vroundpd       $3, %ymm10, %ymm10      
   vroundpd       $3, %ymm11, %ymm11      

   // q * p, high 53 bits
   vmulpd         0x00(ARG1), %ymm8 , %ymm4   
   vmulpd         0x20(ARG1), %ymm9 , %ymm5   
   vmulpd         0x40(ARG1), %ymm10, %ymm6   
   vmulpd         0x60(ARG1), %ymm11, %ymm7   
   
   // q * p, low 53 bits
   vfmsub132pd    0x00(ARG1), %ymm4, %ymm8 
   vfmsub132pd    0x20(ARG1), %ymm5, %ymm9 
   vfmsub132pd    0x40(ARG1), %ymm6, %ymm10
   vfmsub132pd    0x60(ARG1), %ymm7, %ymm11
   
   // (a*b - q*p).hi53
   vsubpd         %ymm4, %ymm0, %ymm4   
   vsubpd         %ymm5, %ymm1, %ymm5   
   vsubpd         %ymm6, %ymm2, %ymm6   
   vsubpd         %ymm7, %ymm3, %ymm7   
                                 
   // (a*b - q*p).lo53           
   vsubpd         %ymm8 , %ymm12, %ymm8    
   vsubpd         %ymm9 , %ymm13, %ymm9    
   vsubpd         %ymm10, %ymm14, %ymm10   
   vsubpd         %ymm11, %ymm15, %ymm11   
   
   // add those hi and lo results, which equal the remainders, up to a possible small
   // multiple of the modulus which will be removed in a final error-correction
   // step. That needs us to preserve both the hi&lo inputs here (since if the
   // sum > 53 bits, e.g. for 52/53-bit inputs and a quotient off by > +-1) we may
   // lose >= 1 one-bits from the bottom):
   vaddpd         %ymm8 , %ymm4, %ymm12   
   vaddpd         %ymm9 , %ymm5, %ymm13   
   vaddpd         %ymm10, %ymm6, %ymm14   
   vaddpd         %ymm11, %ymm7, %ymm15   
   
   // floor((a*b - q*p)/p) gives needed additional multiple of p, e*p, which must
   // be subtracted to correct for the approximateness of q
   vmulpd         0x00(ARG2), %ymm12, %ymm12
   vmulpd         0x20(ARG2), %ymm13, %ymm13
   vmulpd         0x40(ARG2), %ymm14, %ymm14
   vmulpd         0x60(ARG2), %ymm15, %ymm15
   
   vroundpd       $1, %ymm12, %ymm12        
   vroundpd       $1, %ymm13, %ymm13        
   vroundpd       $1, %ymm14, %ymm14        
   vroundpd       $1, %ymm15, %ymm15        
   
   // e * p, high 53 bit
   vmulpd         0x00(ARG1), %ymm12, %ymm0
   vmulpd         0x20(ARG1), %ymm13, %ymm1
   vmulpd         0x40(ARG1), %ymm14, %ymm2
   vmulpd         0x60(ARG1), %ymm15, %ymm3
   
   // e * p, low 53 bits
   vfmsub132pd    0x00(ARG1), %ymm0, %ymm12  
   vfmsub132pd    0x20(ARG1), %ymm1, %ymm13  
   vfmsub132pd    0x40(ARG1), %ymm2, %ymm14  
   vfmsub132pd    0x60(ARG1), %ymm3, %ymm15  
   
   // Error-corrected (a*b - q*p).hi53
   vsubpd         %ymm0, %ymm4, %ymm4   
   vsubpd         %ymm1, %ymm5, %ymm5   
   vsubpd         %ymm2, %ymm6, %ymm6   
   vsubpd         %ymm3, %ymm7, %ymm7   
   
   // Error-corrected (a*b - q*p).lo53
   vsubpd         %ymm12, %ymm8 , %ymm8    
   vsubpd         %ymm13, %ymm9 , %ymm9    
   vsubpd         %ymm14, %ymm10, %ymm10   
   vsubpd         %ymm15, %ymm11, %ymm11   
   
   // Error-corrected remainders
   vaddpd         %ymm8 , %ymm4, %ymm12
   vaddpd         %ymm9 , %ymm5, %ymm13   
   vaddpd         %ymm10, %ymm6, %ymm14  
   vaddpd         %ymm11, %ymm7, %ymm15  
   
   // store a back to memory
   //vmovapd         %ymm12, 0x00(ARG1)        
   //vmovapd         %ymm13, 0x20(ARG1)        
   //vmovapd         %ymm14, 0x40(ARG1)        
   //vmovapd         %ymm15, 0x60(ARG1)
   
   ret
