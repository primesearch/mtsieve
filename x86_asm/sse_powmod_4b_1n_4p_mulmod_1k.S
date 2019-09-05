/* sse_powmod_4b_1n_4p_mulmod_1k.S -- (C) Mark Rodenkirch, January 2018

   void   sse_powmod_4b_1n_4p_mul_1k(uint64_t *b, uint64_t n, uint64_t *p, uint64_t mul);
      Computes b = mul*b^n mod p.
         b is a array of 4 bases.
         p is a array of 4 primes.

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
#define a1     40(%rsp)
#define a2     48(%rsp)
#define a3     56(%rsp)
#define a4     64(%rsp)
#define mul    72(%rsp)
#define old_mxcsr  80(%rsp)
#define new_mxcsr  84(%rsp)
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define ARG4      %rcx
#define t1     -8(%rsp)
#define t2    -16(%rsp)
#define t3    -24(%rsp)
#define t4    -32(%rsp)
#define a1    -40(%rsp)
#define a2    -48(%rsp)
#define a3    -56(%rsp)
#define a4    -64(%rsp)
#define mul   -72(%rsp)
#define old_mxcsr -80(%rsp)
#define new_mxcsr -84(%rsp)
#endif

#define b1      0(%rbp)
#define b2      8(%rbp)
#define b3     16(%rbp)
#define b4     24(%rbp)

#define p1      0(%rdx)
#define p2      8(%rdx)
#define p3     16(%rdx)
#define p4     24(%rdx)

  .globl _sse_powmod_4b_1n_4p_mulmod_1k
  .globl sse_powmod_4b_1n_4p_mulmod_1k

  .p2align 4,,15

_sse_powmod_4b_1n_4p_mulmod_1k:
sse_powmod_4b_1n_4p_mulmod_1k:
    
   push    %rbp
   push    %rbx
   push    %r12
   push    %r13
   push    %r14
   push    %r15
  
#ifdef _WIN64
   push    %rsi
   push    %rdi
   
   sub     $120, %rsp
#endif

   mov     ARG1, %rbp
   mov     ARG2, %rbx
   mov     ARG3, %rdx
   mov     ARG4, mul
   
   stmxcsr old_mxcsr
   mov     old_mxcsr, %eax
   or      $0x6000, %eax      /* Round to zero */
   mov     %eax, new_mxcsr
   ldmxcsr new_mxcsr

   mov     $1, %eax

   cvtsi2sd %eax, %xmm5
   cvtsi2sd %eax, %xmm6
   cvtsi2sd %eax, %xmm7
   cvtsi2sd %eax, %xmm8

   mov     p1, %r8
   mov     p2, %r9
   mov     p3, %r10
   mov     p4, %r11
   
   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4

   // Now we have 1/p for the 4 input values
   divsd   %xmm1, %xmm5
   divsd   %xmm2, %xmm6
   divsd   %xmm3, %xmm7
   divsd   %xmm4, %xmm8
   
   mov     b1, %r8
   mov     b2, %r9
   mov     b3, %r10
   mov     b4, %r11
   
   mov     %r8, a1
   mov     %r9, a2
   mov     %r10, a3
   mov     %r11, a4
   
   cvtsi2sd %r8, %xmm11
   cvtsi2sd %r9, %xmm12
   cvtsi2sd %r10, %xmm13
   cvtsi2sd %r11, %xmm14
   
   bsr     %rbx, %rcx
   dec     %ecx        /* SF=0 predicted */
   jl      vec_done
   mov     $1, %eax
   shl     %cl, %rax   /* second highest bit of n */

/* left-right-powmod(b,n,p)
     a <-- b
     x <-- most significant bit of n
     while x > 0
       b <-- b^2 (mod p)
       x <-- x-1
       if bit x of n is set
         b <-- a*b (mod p)
     return b
*/

  .p2align 4,,1

vec_loop:
   // %rbp -> b1, b2, b3, b4
   // %rsp -> a1, a2, a3, a4
   // %rbx -> n
   // %rax -> current bit of n
   // %rdx -> p1, p2, p3, p4
   // %xmm5-%xmm8 -> 1/p1, 1/p2, 1/p3, 1/p4
   // %xmm11-%xmm14 -> a1, a2, a3, a4

  .p2align 4,,7

vec_sqr:
   // Compute b^2 mod p for each b
   mov     b1, %r8
   mov     b2, %r9
   mov     b3, %r10
   mov     b4, %r11
   
   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4

   mulsd   %xmm1, %xmm1
   mulsd   %xmm2, %xmm2
   mulsd   %xmm3, %xmm3
   mulsd   %xmm4, %xmm4
   
   mulsd   %xmm5, %xmm1
   mulsd   %xmm6, %xmm2
   mulsd   %xmm7, %xmm3
   mulsd   %xmm8, %xmm4

   cvtsd2si %xmm1, %r12       // r12 = trunc(b1^2 / p1)
   cvtsd2si %xmm2, %r13       // r13 = trunc(b2^2 / p2)
   cvtsd2si %xmm3, %r14       // r14 = trunc(b3^2 / p3)
   cvtsd2si %xmm4, %r15       // r15 = trunc(b4^2 / p4)

   // Yes, this can be a 104-bit product, but we only care
   // about the lower 52 bits because p < 2^52.
   imul    %r8, %r8
   imul    %r9, %r9
   imul    %r10, %r10
   imul    %r11, %r11

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
   
   mov     %r8, b1            // b1 = b1*b1 mod p1
   mov     %r9, b2            // b2 = b2*b2 mod p2
   mov     %r10, b3           // b3 = b3*b3 mod p3
   mov     %r11, b4           // b4 = b4*b4 mod p4

   sub     p1, %r8            // if b1 > p1, then b1 -= p1
   jc      s0
   mov     %r8, b1
   
s0:
   sub     p2, %r9            // if b2 > p2, then b2 -= p2
   jc      s1
   mov     %r9, b2
   
s1:
   sub     p3, %r10           // if b3 > p3, then b3 -= p3
   jc      s2
   mov     %r10, b3
   
s2:
   sub     p4, %r11           // if b4 > p4, then b4 -= p4
   jc      s3
   mov     %r11, b4
   
s3:
   test    %rax, %rbx         /* ZF unpredictable */
   jz      vec_shift

  .p2align 4,,7

vec_mul:   
   // Compute a*b mod p for each b
   mov     b1, %r8
   mov     b2, %r9
   mov     b3, %r10
   mov     b4, %r11
   
   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4
   
   mulsd   %xmm11, %xmm1
   mulsd   %xmm12, %xmm2
   mulsd   %xmm13, %xmm3
   mulsd   %xmm14, %xmm4

   mulsd   %xmm5, %xmm1
   mulsd   %xmm6, %xmm2
   mulsd   %xmm7, %xmm3
   mulsd   %xmm8, %xmm4
   
   cvtsd2si %xmm1, %r12       // r12 = trunc(b1^a1 / p1)
   cvtsd2si %xmm2, %r13       // r13 = trunc(b2^a2 / p2)
   cvtsd2si %xmm3, %r14       // r14 = trunc(b3^a3 / p3)
   cvtsd2si %xmm4, %r15       // r15 = trunc(b4^a4 / p4)
    
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
   
   mov     %r8, b1            // b1 = a1*b1 mod p1
   mov     %r9, b2            // b2 = a2*b2 mod p2
   mov     %r10, b3           // b3 = a3*b3 mod p3
   mov     %r11, b4           // b4 = a4*b4 mod p4

   sub     p1, %r8            // if b1 > p1, then b1 -= p1
   jc      m0
   mov     %r8, b1
   
m0:
   sub     p2, %r9            // if b2 > p2, then b2 -= p2
   jc      m1
   mov     %r9, b2
   
m1:
   sub     p3, %r10           // if b3 > p3, then b3 -= p3
   jc      m2
   mov     %r10, b3
   
m2:
   sub     p4, %r11           // if b4 > p4, then b4 -= p4
   jc      m3
   mov     %r11, b4
   
m3:
   .p2align 4,,2

vec_shift:
   shr     %rax
   jnz     vec_loop

vec_done:
   // we have b^n (mod p), now compute mul*b^n (mod p)
   mov     b1, %r8
   mov     b2, %r9
   mov     b3, %r10
   mov     b4, %r11
   
   cvtsi2sd %r8, %xmm1
   cvtsi2sd %r9, %xmm2
   cvtsi2sd %r10, %xmm3
   cvtsi2sd %r11, %xmm4
   
   mov     mul, %r8
   mov     mul, %r9
   mov     mul, %r10
   mov     mul, %r11
   
   cvtsi2sd %r8, %xmm11
   cvtsi2sd %r9, %xmm12
   cvtsi2sd %r10, %xmm13
   cvtsi2sd %r11, %xmm14

   mulsd   %xmm11, %xmm1
   mulsd   %xmm12, %xmm2
   mulsd   %xmm13, %xmm3
   mulsd   %xmm14, %xmm4

   mulsd   %xmm5, %xmm1
   mulsd   %xmm6, %xmm2
   mulsd   %xmm7, %xmm3
   mulsd   %xmm8, %xmm4
   
   cvtsd2si %xmm1, %r12       // r12 = trunc(b1^a1 / p1)
   cvtsd2si %xmm2, %r13       // r13 = trunc(b2^a2 / p2)
   cvtsd2si %xmm3, %r14       // r14 = trunc(b3^a3 / p3)
   cvtsd2si %xmm4, %r15       // r15 = trunc(b4^a4 / p4)
    
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
   
   mov     %r8, b1            // b1 = a1*b1 mod p1
   mov     %r9, b2            // b2 = a2*b2 mod p2
   mov     %r10, b3           // b3 = a3*b3 mod p3
   mov     %r11, b4           // b4 = a4*b4 mod p4

   sub     p1, %r8            // if b1 > p1, then b1 -= p1
   jc      lm0
   mov     %r8, b1
   
lm0:
   sub     p2, %r9            // if b2 > p2, then b2 -= p2
   jc      lm1
   mov     %r9, b2
   
lm1:
   sub     p3, %r10           // if b3 > p3, then b3 -= p3
   jc      lm2
   mov     %r10, b3
   
lm2:
   sub     p4, %r11           // if b4 > p4, then b4 -= p4
   jc      lm3
   mov     %r11, b4
   
lm3:
   ldmxcsr old_mxcsr

#ifdef _WIN64
   add     $120, %rsp
   
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
