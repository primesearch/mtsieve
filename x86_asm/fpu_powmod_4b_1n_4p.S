/* fpu_powmod_4b_1n_4p.S -- (C) Mark Rodenkirch, January 2018

   void   fpu_powmod_4b_1n_4p(uint64_t *b, uint64_t n, uint64_t *p);
      Computes b = b^n mod p.
         b is a array of 4 bases.
         p is a array of 4 primes.

      Assumes FPU is set to double extended precision and round to zero.

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
#define a1     40(%rsp)
#define a2     48(%rsp)
#define a3     56(%rsp)
#define a4     64(%rsp)
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define t1     -8(%rsp)
#define t2    -16(%rsp)
#define t3    -24(%rsp)
#define t4    -32(%rsp)
#define a1    -40(%rsp)
#define a2    -48(%rsp)
#define a3    -56(%rsp)
#define a4    -64(%rsp)
#endif

#define b1      0(%rbp)
#define b2      8(%rbp)
#define b3     16(%rbp)
#define b4     24(%rbp)

#define p1      0(%rdx)
#define p2      8(%rdx)
#define p3     16(%rdx)
#define p4     24(%rdx)

  .globl _fpu_powmod_4b_1n_4p
  .globl fpu_powmod_4b_1n_4p

  .p2align 4,,15

_fpu_powmod_4b_1n_4p:
fpu_powmod_4b_1n_4p:
    
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
   
// Compute 1.0 / p for each p
   fildll  p4
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p4
   
   fildll  p3
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p3
   
   fildll  p2
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p2
   
   fildll  p1
   fld1                       // put 1.0 into %st(0)
   fdivp   %st(0), %st(1)     // compute 1.0 / p1
   
   mov     b1, %r8
   mov     b2, %r9
   mov     b3, %r10
   mov     b4, %r11
      
   mov     %r8, a1
   mov     %r9, a2
   mov     %r10, a3
   mov     %r11, a4
      
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
   // %st(0,1,2,3) -> 1/p1, 1/p2, 1/p3, 1/p4

  .p2align 4,,7

vec_sqr:
   // Compute b^2 mod p for each b
   fildll  b1
   fmul    %st(0), %st(0)
   fmul    %st(1), %st(0)
   fistpll t1                 // t1 = trunc(b1^2 / p1)

   fildll  b2
   fmul    %st(0), %st(0)
   fmul    %st(2), %st(0)
   fistpll t2                 // t2 = trunc(b2^2 / p2)

   fildll  b3
   fmul    %st(0), %st(0)
   fmul    %st(3), %st(0)
   fistpll t3                 // t3 = trunc(b3^2 / p3)
   
   fildll  b4
   fmul    %st(0), %st(0)
   fmul    %st(4), %st(0)
   fistpll t4                 // t4 = trunc(b4^2 / p4)

   mov     b1, %r8
   mov     b2, %r9
   mov     b3, %r10
   mov     b4, %r11
   
   // Yes, this can be a 124-bit product, but we only care
   // about the lower 62 bits because p < 2^62.
   imul    %r8, %r8
   imul    %r9, %r9
   imul    %r10, %r10
   imul    %r11, %r11

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
   
   mov     %r8, b1            // b1 = b1*b1 mod p1
   mov     %r9, b2            // b2 = b2*b2 mod p2
   mov     %r10, b3           // b3 = b3*b3 mod p3
   mov     %r11, b4           // b4 = b4*b4 mod p4

   sub     p1, %r8            // if b1 > p1, then b1 -= p1
   jl      s0
   mov     %r8, b1
   
s0:
   sub     p2, %r9            // if b2 > p2, then b2 -= p2
   jl      s1
   mov     %r9, b2
   
s1:
   sub     p3, %r10           // if b3 > p3, then b3 -= p3
   jl      s2
   mov     %r10, b3
   
s2:
   sub     p4, %r11           // if b4 > p4, then b4 -= p4
   jl      s3
   mov     %r11, b4
   
s3:
   test    %rax, %rbx         /* ZF unpredictable */
   jz      vec_shift

  .p2align 4,,7

vec_mul:
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
   
   mov     %r8, b1            // b1 = a1*b1 mod p1
   mov     %r9, b2            // b2 = a2*b2 mod p2
   mov     %r10, b3           // b3 = a3*b3 mod p3
   mov     %r11, b4           // b4 = a4*b4 mod p4

   sub     p1, %r8            // if b1 > p1, then b1 -= p1
   jl      m0
   mov     %r8, b1
   
m0:
   sub     p2, %r9            // if b2 > p2, then b2 -= p2
   jl      m1
   mov     %r9, b2
   
m1:
   sub     p3, %r10           // if b3 > p3, then b3 -= p3
   jl      m2
   mov     %r10, b3
   
m2:
   sub     p4, %r11           // if b4 > p4, then b4 -= p4
   jl      m3
   mov     %r11, b4
   
m3:
   .p2align 4,,2

vec_shift:
   shr     %rax
   jnz     vec_loop

vec_done:
   fstp    %st(0)
   fstp    %st(0)
   fstp    %st(0)
   fstp    %st(0)

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
