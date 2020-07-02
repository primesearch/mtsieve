/* avx_powmod.S -- (C) Mark Rodenkirch, May 2018

   All AVX functions require parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[16];

   uint16_t  avx_powmod(double *b, uint64_t n, double *p, double *reciprocal)
      Compute b = b^n mod p for 1 < b,p < 2^52 and 1 < n < 2^64
      
      b, p, and reciprocal are arrays of 16 doubles.
      
      You must call avx_mod_get_inverse to compute the reciprocal that is passed
      to this function.

      Upon output b will contain the 16 remainders of the powmod.

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
#define ARG4      %r9

#define a1   0x00(%rsp)
#define a2   0x20(%rsp)
#define a3   0x40(%rsp)
#define a4   0x60(%rsp)
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define ARG4      %rcx

#define a1  -0x00(%rsp)
#define a2  -0x20(%rsp)
#define a3  -0x40(%rsp)
#define a4  -0x60(%rsp)
#endif

#define b1   0x00(%r10)
#define b2   0x20(%r10)
#define b3   0x40(%r10)
#define b4   0x60(%r10)

#define EXP       %r9

#define p1   0x00(%r8)
#define p2   0x20(%r8)
#define p3   0x40(%r8)
#define p4   0x60(%r8)

#define inv1 0x00(%r11)
#define inv2 0x20(%r11)
#define inv3 0x40(%r11)
#define inv4 0x60(%r11)

  .p2align 4,,15
  
  .globl _avx_powmod
  .globl avx_powmod

 _avx_powmod:
 avx_powmod:

 // We need %rcx below, so copy it to another register
   // that does not require a push/pop
   mov            ARG1, %r10
   mov            ARG4, %r11
   mov            ARG2, EXP
   mov            ARG3, %r8

   // Ensure %rsp is on a 32-byte boundary
   mov            %rsp, %rdx

#ifdef _WIN64
   sub            $400, %rsp
rsp_loop:
   test           $31, %rsp
   jz             rsp_done
   sub            $1, %rsp
   jmp            rsp_loop
#else
   sub            $400, %rsp
rsp_loop:
   test           $31, %rsp
   jz             rsp_done
   add            $1, %rsp
   jmp            rsp_loop
#endif
   
rsp_done:
   vmovapd        b1, %ymm12
   vmovapd        b2, %ymm13
   vmovapd        b3, %ymm14
   vmovapd        b4, %ymm15

   // If n = 1, then we are done
   bsr     EXP, %rcx
   dec     %rcx
   jl      all_done

   // Note that %cl is the lowest 8 bits of %rcx
   mov     $1, %rax
   shl     %cl, %rax

   // left-right-powmod(b,n,p)
   //   a <-- b
   //   x <-- most significant bit of n
   //   while x > 0
   //     b <-- b^2 (mod p)
   //     x <-- x-1
   //     if bit x of n is set
   //       b <-- a*b (mod p)
   //   return b

   // Fingers crossed that %rsp is aligned to 32-byte boundary
   // otherwise have to write code to ensure that happens.
   vmovapd        %ymm12, a1
   vmovapd        %ymm13, a2
   vmovapd        %ymm14, a3
   vmovapd        %ymm15, a4
   
sqr_step:
   vmovapd        %ymm12, %ymm0
   vmovapd        %ymm13, %ymm1
   vmovapd        %ymm14, %ymm2
   vmovapd        %ymm15, %ymm3

   // b * b, high 53 bits
   vmulpd		   %ymm12, %ymm0, %ymm4
   vmulpd		   %ymm13, %ymm1, %ymm5
   vmulpd		   %ymm14, %ymm2, %ymm6
   vmulpd		   %ymm15, %ymm3, %ymm7
   
   // q = trunc(b * b / p)
   vmulpd		   inv1, %ymm4, %ymm8    
   vmulpd		   inv2, %ymm5, %ymm9    
   vmulpd		   inv3, %ymm6, %ymm10   
   vmulpd		   inv4, %ymm7, %ymm11   
   
   // b * b, low 53 bits
   vfmsub132pd    %ymm0, %ymm4, %ymm12  
   vfmsub132pd    %ymm1, %ymm5, %ymm13  
   vfmsub132pd    %ymm2, %ymm6, %ymm14  
   vfmsub132pd    %ymm3, %ymm7, %ymm15  
   
   // Round q toward 0 (since these are >= 0 could also round toward -oo)
   vroundpd	      $3, %ymm8 , %ymm8       
   vroundpd	      $3, %ymm9 , %ymm9       
   vroundpd	      $3, %ymm10, %ymm10      
   vroundpd	      $3, %ymm11, %ymm11      

   // q * p, high 53 bits
   vmulpd		   p1, %ymm8 , %ymm0   
   vmulpd		   p2, %ymm9 , %ymm1   
   vmulpd		   p3, %ymm10, %ymm2   
   vmulpd		   p4, %ymm11, %ymm3   
   
   // q * p, low 53 bits
   vfmsub132pd	   p1, %ymm0, %ymm8 
   vfmsub132pd	   p2, %ymm1, %ymm9 
   vfmsub132pd	   p3, %ymm2, %ymm10
   vfmsub132pd	   p4, %ymm3, %ymm11
   
   // (b*b - q*p).hi53
   vsubpd		   %ymm0, %ymm4, %ymm0   
   vsubpd		   %ymm1, %ymm5, %ymm1   
   vsubpd		   %ymm2, %ymm6, %ymm2   
   vsubpd		   %ymm3, %ymm7, %ymm3   
                                 
   // (b*b - q*p).lo53           
   vsubpd		   %ymm8 , %ymm12, %ymm8    
   vsubpd		   %ymm9 , %ymm13, %ymm9    
   vsubpd		   %ymm10, %ymm14, %ymm10   
   vsubpd		   %ymm11, %ymm15, %ymm11   
   
   // add those hi and lo results, which equal the remainders, up to a possible small
   // multiple of the modulus which will be removed in a final error-correction
	// step. That needs us to preserve both the hi&lo inputs here (since if the
   // sum > 53 bits, e.g. for 52/53-bit inputs and a quotient off by > +-1) we may
   // lose >= 1 one-bits from the bottom):
   vaddpd		   %ymm8 , %ymm0, %ymm12   
   vaddpd		   %ymm9 , %ymm1, %ymm13   
   vaddpd		   %ymm10, %ymm2, %ymm14   
   vaddpd		   %ymm11, %ymm3, %ymm15   
   
   // floor((b*b - q*p)/p) gives needed additional multiple of p, e*p, which must
   // be subtracted to correct for the approximateness of q
   vmulpd		   inv1, %ymm12, %ymm12
   vmulpd		   inv2, %ymm13, %ymm13
   vmulpd		   inv3, %ymm14, %ymm14
   vmulpd		   inv4, %ymm15, %ymm15
   
   vroundpd       $1, %ymm12, %ymm12        
   vroundpd       $1, %ymm13, %ymm13        
   vroundpd       $1, %ymm14, %ymm14        
   vroundpd       $1, %ymm15, %ymm15        
   
   // e * p, high 53 bit
   vmulpd		   p1, %ymm12, %ymm4
   vmulpd		   p2, %ymm13, %ymm5
   vmulpd		   p3, %ymm14, %ymm6
   vmulpd		   p4, %ymm15, %ymm7
   
   // e * p, low 53 bits
   vfmsub132pd	   p1, %ymm4, %ymm12  
   vfmsub132pd	   p2, %ymm5, %ymm13  
   vfmsub132pd	   p3, %ymm6, %ymm14  
   vfmsub132pd	   p4, %ymm7, %ymm15  
   
   // Error-corrected (b*b - q*p).hi53
   vsubpd		   %ymm4, %ymm0, %ymm0   
   vsubpd		   %ymm5, %ymm1, %ymm1   
   vsubpd		   %ymm6, %ymm2, %ymm2   
   vsubpd		   %ymm7, %ymm3, %ymm3   
   
   // Error-corrected (b*b - q*p).lo53
   vsubpd		   %ymm12, %ymm8 , %ymm8    
   vsubpd		   %ymm13, %ymm9 , %ymm9    
   vsubpd		   %ymm14, %ymm10, %ymm10   
   vsubpd		   %ymm15, %ymm11, %ymm11   
   
   // Error-corrected remainders
   vaddpd		   %ymm8 , %ymm0, %ymm12
   vaddpd		   %ymm9 , %ymm1, %ymm13   
   vaddpd		   %ymm10, %ymm2, %ymm14
   vaddpd		   %ymm11, %ymm3, %ymm15  

   // ymm12:ymm15 contain b where b = b^2 mod p
   
   // If the current bit is 0, then skip the multiply step
   test           %rax, EXP
   jz             bit_shift
   
mul_step:
   // a * b, high 53 bits
   vmulpd		   a1, %ymm12, %ymm4
   vmulpd		   a2, %ymm13, %ymm5
   vmulpd		   a3, %ymm14, %ymm6
   vmulpd		   a4, %ymm15, %ymm7
   
   // q = trunc(a * b / p)
   vmulpd		   inv1, %ymm4, %ymm8    
   vmulpd		   inv2, %ymm5, %ymm9    
   vmulpd		   inv3, %ymm6, %ymm10   
   vmulpd		   inv4, %ymm7, %ymm11   
   
   // a * b, low 53 bits
   vfmsub132pd    a1, %ymm4, %ymm12  
   vfmsub132pd    a2, %ymm5, %ymm13  
   vfmsub132pd    a3, %ymm6, %ymm14  
   vfmsub132pd    a4, %ymm7, %ymm15  
   
   // Round q toward 0 (since these are >= 0 could also round toward -oo)
   vroundpd	      $3, %ymm8 , %ymm8       
   vroundpd	      $3, %ymm9 , %ymm9       
   vroundpd	      $3, %ymm10, %ymm10      
   vroundpd	      $3, %ymm11, %ymm11      

   // q * p, high 53 bits
   vmulpd		   p1, %ymm8 , %ymm0   
   vmulpd		   p2, %ymm9 , %ymm1   
   vmulpd		   p3, %ymm10, %ymm2   
   vmulpd		   p4, %ymm11, %ymm3   
   
   // q * p, low 53 bits
   vfmsub132pd	   p1, %ymm0, %ymm8 
   vfmsub132pd	   p2, %ymm1, %ymm9 
   vfmsub132pd	   p3, %ymm2, %ymm10
   vfmsub132pd	   p4, %ymm3, %ymm11
   
   // (a*b - q*p).hi53
   vsubpd		   %ymm0, %ymm4, %ymm0   
   vsubpd		   %ymm1, %ymm5, %ymm1   
   vsubpd		   %ymm2, %ymm6, %ymm2   
   vsubpd		   %ymm3, %ymm7, %ymm3   
                                 
   // (a*b - q*p).lo53           
   vsubpd		   %ymm8 , %ymm12, %ymm8    
   vsubpd		   %ymm9 , %ymm13, %ymm9    
   vsubpd		   %ymm10, %ymm14, %ymm10   
   vsubpd		   %ymm11, %ymm15, %ymm11   
   
   // add those hi and lo results, which equal the remainders, up to a possible small
   // multiple of the modulus which will be removed in a final error-correction
	// step. That needs us to preserve both the hi&lo inputs here (since if the
   // sum > 53 bits, e.g. for 52/53-bit inputs and a quotient off by > +-1) we may
   // lose >= 1 one-bits from the bottom):
   vaddpd		   %ymm8 , %ymm0, %ymm12   
   vaddpd		   %ymm9 , %ymm1, %ymm13   
   vaddpd		   %ymm10, %ymm2, %ymm14   
   vaddpd		   %ymm11, %ymm3, %ymm15   
   
   // floor((a*b - q*p)/p) gives needed additional multiple of p, e*p, which must
   // be subtracted to correct for the approximateness of q
   vmulpd		   inv1, %ymm12, %ymm12
   vmulpd		   inv2, %ymm13, %ymm13
   vmulpd		   inv3, %ymm14, %ymm14
   vmulpd		   inv4, %ymm15, %ymm15
   
   vroundpd       $1, %ymm12, %ymm12        
   vroundpd       $1, %ymm13, %ymm13        
   vroundpd       $1, %ymm14, %ymm14        
   vroundpd       $1, %ymm15, %ymm15        
   
   // e * p, high 53 bit
   vmulpd		   p1, %ymm12, %ymm4
   vmulpd		   p2, %ymm13, %ymm5
   vmulpd		   p3, %ymm14, %ymm6
   vmulpd		   p4, %ymm15, %ymm7
   
   // e * p, low 53 bits
   vfmsub132pd	   p1, %ymm4, %ymm12  
   vfmsub132pd	   p2, %ymm5, %ymm13  
   vfmsub132pd	   p3, %ymm6, %ymm14  
   vfmsub132pd	   p4, %ymm7, %ymm15  
   
   // Error-corrected (a*b - q*p).hi53
   vsubpd		   %ymm4, %ymm0, %ymm0   
   vsubpd		   %ymm5, %ymm1, %ymm1   
   vsubpd		   %ymm6, %ymm2, %ymm2   
   vsubpd		   %ymm7, %ymm3, %ymm3   
   
   // Error-corrected (a*b - q*p).lo53
   vsubpd		   %ymm12, %ymm8 , %ymm8    
   vsubpd		   %ymm13, %ymm9 , %ymm9    
   vsubpd		   %ymm14, %ymm10, %ymm10   
   vsubpd		   %ymm15, %ymm11, %ymm11   
   
   // Error-corrected remainders
   vaddpd		   %ymm8 , %ymm0, %ymm12
   vaddpd		   %ymm9 , %ymm1, %ymm13   
   vaddpd		   %ymm10, %ymm2, %ymm14  
   vaddpd		   %ymm11, %ymm3, %ymm15  
   
   // ymm12:ymm15 contain b where b = a*b mod p
   
bit_shift:
   shr            %rax
   jnz            sqr_step
   
all_done:
   vmovapd        %ymm12, b1
   vmovapd        %ymm13, b2
   vmovapd        %ymm14, b3
   vmovapd        %ymm15, b4
  
   mov            %rdx, %rsp
   
   ret
   
