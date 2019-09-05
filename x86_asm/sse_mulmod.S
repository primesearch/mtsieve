/* sse_mulmod.S -- (C) Mark Rodenkirch, May 2018

   void   sse_mulmod(uint64_t a, uint64_t b, uint64_t p, double *bdivp);
      Computes a = a*b mod p.
         bdivp points to b/p

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
#else
#define ARG1      %rdi
#define ARG2      %rsi
#define ARG3      %rdx
#define ARG4      %rcx
#endif

  .globl _sse_mulmod
  .globl sse_mulmod

  .p2align 4,,15

_sse_mulmod:
sse_mulmod:

   push    %r12

   cvtsi2sd ARG1, %xmm1
   
   movsd   0(ARG4), %xmm5
   
   // Compute a * b/p
   mulsd   %xmm5, %xmm1
   
   cvtsd2si %xmm1, %r12       // r12 = trunc(a1*b1 / p1)
   
   // Yes, this can be a 104-bit product, but we only care
   // about the lower 52 bits because p < 2^52.
   imul    ARG1, ARG2
   
   // Yes, this can be a 104-bit product, but we only care
   // about the lower 52 bits because p < 2^52.
   imul    ARG3, %r12
   
   // Since p < 2^52, we know that the upper 52 bits of each
   // product is the same, so they have no impact on the result
   sub     %r12, ARG2
   
   mov     ARG2, %rax         // a1 = a1*b1 mod p1
   sub     ARG3, ARG2         // if a > p, then a -= p
   jc      m0
   mov     ARG2, %rax
   
m0:
   .p2align 4,,2
   
   pop     %r12
   
   ret
