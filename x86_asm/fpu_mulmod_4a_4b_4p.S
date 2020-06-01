/* fpu_mulmod_4a_4b_4p.S -- (C) Mark Rodenkirch, May 2018

   void   fpu_mulmod_4a_4b_4p(uint64_t *a, uint64_t *b, uint64_t *p);
      Computes a = a*b mod p.
         a is an array of 4 numbers.
         b is an array of 4 numbers.
         p is an array of 4 primes.

      Assumes FPU is set to double extended precision and round to zero

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
#define t1      8(%rsp)
#define t2     16(%rsp)
#define t3     24(%rsp)
#define t4     32(%rsp)
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define t1     -8(%rsp)
#define t2    -16(%rsp)
#define t3    -24(%rsp)
#define t4    -32(%rsp)
#endif

#define a1      0(%rbp)
#define a2      8(%rbp)
#define a3     16(%rbp)
#define a4     24(%rbp)

#define b1      0(%rbx)
#define b2      8(%rbx)
#define b3     16(%rbx)
#define b4     24(%rbx)

#define p1      0(%rdx)
#define p2      8(%rdx)
#define p3     16(%rdx)
#define p4     24(%rdx)

  .globl _fpu_mulmod_4a_4b_4p
  .globl fpu_mulmod_4a_4b_4p

  .p2align 4,,15

_fpu_mulmod_4a_4b_4p:
fpu_mulmod_4a_4b_4p:
    
   push    %rbp
   push    %rbx
   push    %r12
   push    %r13
   push    %r14
   push    %r15
  
#ifdef _WIN64
   push    %rsi
   push    %rdi
   
   sub     $60, %rsp
#endif

   mov     ARG1, %rbp
   mov     ARG2, %rbx
   mov     ARG3, %rdx

   // Compute a*b mod p for each b
   fildll  a1
   fildll  b1
   fmulp   %st(0), %st(1)
   fmul    %st(1), %st(0)
   fistpll t1                 // t1 = trunc(a1*b1 / p1)

   fildll  a2
   fildll  b2
   fmulp   %st(0), %st(1)
   fmul    %st(2), %st(0)
   fistpll t2                 // t2 = trunc(a2*b2 / p2)

   fildll  a3
   fildll  b3
   fmulp   %st(0), %st(1)
   fmul    %st(3), %st(0)
   fistpll t3                 // t3 = trunc(a3*b3 / p3)
   
   fildll  a4
   fildll  b4
   fmulp   %st(0), %st(1)
   fmul    %st(4), %st(0)
   fistpll t4                 // t4 = trunc(a4*b4 / p4)
      
   mov     a1, %r8
   mov     a2, %r9
   mov     a3, %r10
   mov     a4, %r11
   
   // Yes, this can be a 124-bit product, but we only care
   // about the lower 62 bits because p < 2^62.
   imul    b1, %r8
   imul    b2, %r9
   imul    b3, %r10
   imul    b4, %r11

   mov     t1, %r12
   mov     t2, %r13
   mov     t3, %r14
   mov     t4, %r15
   
   // Yes, this can be a 124-bit product, but we only care
   // about the lower 62 bits because p < 2^62.
   imul    p1, %r12
   imul    p2, %r13
   imul    p3, %r14
   imul    p4, %r15
   
   // Since p < 2^62, we know that the upper 62 bits of each
   // product is the same, so they have no impact on the result
   sub     %r12, %r8
   sub     %r13, %r9
   sub     %r14, %r10
   sub     %r15, %r11
   
   mov     %r8, a1            // a1 = a1*b1 mod p1
   mov     %r9, a2            // a2 = a2*b2 mod p2
   mov     %r10, a3           // a3 = a3*b3 mod p3
   mov     %r11, a4           // a4 = a4*b4 mod p4

   sub     p1, %r8            // if a1 > p1, then a1 -= p1
   jl      m0
   mov     %r8, a1
   
m0:
   sub     p2, %r9            // if a2 > p2, then a2 -= p2
   jl      m1
   mov     %r9, a2
   
m1:
   sub     p3, %r10           // if a3 > p3, then a3 -= p3
   jl      m2
   mov     %r10, a3
   
m2:
   sub     p4, %r11           // if a4 > p4, then a4 -= p4
   jl      m3
   mov     %r11, a4
   
m3:

#ifdef _WIN64
   add     $60, %rsp
   
   pop     %rdi
   pop     %rsi
#endif

   pop     %r15
   pop     %r14
   pop     %r13
   pop     %r12
   pop     %rbx
   pop     %rbp
   
   ret
