/* fpu_mulmod_iter.S -- (C) Mark Rodenkirch, January 2018

   uint64_t fpu_mulmod_iter(uint64_t a, uint64_t b, uint64_t p);
      Returns a*b (mod p), where 0 <= a,b < p < 2^62.

      Assumes FPU is set to double extended precision and round to zero.
      Assumes that st(0) is b/p.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifdef _WIN64
#define ARGa     %rcx
#define ARGb     %rdx
#define ARGp     %r8
#define t1    16(%rsp)
#else
#define ARGa     %rdi
#define ARGb     %rsi
#define ARGp     %rdx
#define t1   -32(%rsp)
#endif

   .text
   .globl _fpu_mulmod_iter
   .globl fpu_mulmod_iter

   .p2align 4,,15
       
_fpu_mulmod_iter:
fpu_mulmod_iter:

#ifdef _WIN64
   # Windows build crashes without this, but fpu_mulmod doesn't need it.  Very odd.
   sub    $40, %rsp
#endif

   push    %r10
   push    %r11

   mov     ARGa, t1
   
   fildll  t1

   fmul    %st(1), %st(0)
   
   fistpll t1                 // t1 = trunc(a*b / p)
   
   mov     ARGa, %r10

   // Yes, this can be a 124-bit product, but we only care
   // about the lower 62 bits because a, b, p < 2^62.
   imul    ARGb, %r10

   mov     t1, %r11
   
   // Yes, this can be a 124-bit product, but we only care
   // about the lower 62 bits because a, b, p < 2^62.
   imul    ARGp, %r11
   
   // Since p < 2^62, we know that the upper 62 bits of each
   // product is the same, so they have no impact on the result
   sub     %r11, %r10
   
   mov     %r10, %rax         // rax = a*b mod p

   sub     ARGp, %r10  /* CF=0 predicted */
   jc      done
   mov     %r10, %rax

done:
   pop     %r11
   pop     %r10
   
#ifdef _WIN64
   add     $40, %rsp
#endif

   ret