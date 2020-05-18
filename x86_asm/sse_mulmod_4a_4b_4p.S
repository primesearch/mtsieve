/* sse_mulmod_4a_4b_4p.S -- (C) Mark Rodenkirch, May 2018

   void   sse_mulmod_4a_4b_4p(uint64_t *a, uint64_t *b, uint64_t *p, double *bdivp);
      Computes a = a*b mod p.
         a is an array of 4 numbers.
         b is an array of 4 numbers.
         p is an array of 4 primes.
         bdivp is an array of b/p

      Assumes SSE is set to round to zero.

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
#define t1      8(%rsp)
#define t2     16(%rsp)
#define t3     24(%rsp)
#define t4     32(%rsp)
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define ARG4      %rcx
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

#define bdivp1  0(%rcx)
#define bdivp2  8(%rcx)
#define bdivp3 16(%rcx)
#define bdivp4 24(%rcx)

  .globl _sse_mulmod_4a_4b_4p
  .globl sse_mulmod_4a_4b_4p

  .p2align 4,,15

_sse_mulmod_4a_4b_4p:
sse_mulmod_4a_4b_4p:
    
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
   mov     ARG4, %rcx

   movsd   bdivp1, %xmm5
   movsd   bdivp2, %xmm6
   movsd   bdivp3, %xmm7
   movsd   bdivp4, %xmm8

   mov     a1, %r8
   mov     a2, %r9
   mov     a3, %r10
   mov     a4, %r11

   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4
      
   // Compute a*b mod p for each b
   mulsd   %xmm5, %xmm1
   mulsd   %xmm6, %xmm2
   mulsd   %xmm7, %xmm3
   mulsd   %xmm8, %xmm4
   
   cvtsd2si %xmm1, %r12       // r12 = trunc(a1*b1 / p1)
   cvtsd2si %xmm2, %r13       // r13 = trunc(a2*b2 / p2)
   cvtsd2si %xmm3, %r14       // r14 = trunc(a3*b3 / p3)
   cvtsd2si %xmm4, %r15       // r15 = trunc(a4*b4 / p4)
   
   mov     a1, %r8
   mov     a2, %r9
   mov     a3, %r10
   mov     a4, %r11
   
   // Yes, this can be a 104-bit product, but we only care
   // about the lower 52 bits because p < 2^52.
   imul    b1, %r8
   imul    b2, %r9
   imul    b3, %r10
   imul    b4, %r11
   
   // Yes, this can be a 104-bit product, but we only care
   // about the lower 52 bits because p < 2^52.
   imul    p1, %r12
   imul    p2, %r13
   imul    p3, %r14
   imul    p4, %r15
   
   // Since p < 2^52, we know that the upper 52 bits of each
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
   .p2align 4,,2
   
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
