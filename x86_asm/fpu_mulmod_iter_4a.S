/* fpu_mulmod_iter_4a.S -- (C) Mark Rodenkirch, January 2018

   uint64_t fpu_mulmod_iter_4a(uint64_t *a, uint64_t b, uint64_t p);
      Returns a*b (mod p), where 0 <= a,b < p < 2^62 for 4 distinct a.

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
#define VAR1  56(%rsp)
#define VAR2  64(%rsp)
#define VAR3  72(%rsp)
#define VAR4  80(%rsp)
#else
#define ARGa     %rdi
#define ARGb     %rsi
#define ARGp     %rdx
#define VAR1 -32(%rsp)
#define VAR2 -24(%rsp)
#define VAR3 -16(%rsp)
#define VAR4  -8(%rsp)
#endif

   .text
   .globl _fpu_mulmod_iter_4a
   .globl fpu_mulmod_iter_4a

   .p2align 4,,15
       
_fpu_mulmod_iter_4a:
fpu_mulmod_iter_4a:

#ifdef _WIN64 
   sub     $120, %rsp
#endif

   push    %rax
   push    %rbx
   push    %r10
   push    %r11
   push    %r12
   push    %r13
   push    %r14
   push    %r15

   mov     0(ARGa), %r10
   mov     8(ARGa), %r11
   mov     16(ARGa), %r12
   mov     24(ARGa), %r13
   
   fildll  0(ARGa)
   fmul    %st(1), %st(0)
   fistpll VAR1
   
   fildll  8(ARGa)
   fmul    %st(1), %st(0)
   fistpll VAR2
   
   fildll  16(ARGa)
   fmul    %st(1), %st(0)
   fistpll VAR3
   
   fildll  24(ARGa)
   fmul    %st(1), %st(0)
   fistpll VAR4

   imul    ARGb, %r10
   imul    ARGb, %r11
   imul    ARGb, %r12
   imul    ARGb, %r13

   mov     VAR1, %rax
   mov     VAR2, %rbx
   mov     VAR3, %r14
   mov     VAR4, %r15
          
   imul    ARGp, %rax
   imul    ARGp, %rbx
   imul    ARGp, %r14
   imul    ARGp, %r15

   sub     %rax, %r10
   sub     %rbx, %r11
   sub     %r14, %r12
   sub     %r15, %r13
   
   mov     %r10, 0(ARGa)
   mov     %r11, 8(ARGa)
   mov     %r12, 16(ARGa)
   mov     %r13, 24(ARGa)
   
   sub     ARGp, %r10
   jl      a1
   mov     %r10, 0(ARGa)
  
a1:
   sub     ARGp, %r11
   jl      a2
   mov     %r11, 8(ARGa)
           
a2:         
   sub     ARGp, %r12
   jl      a3
   mov     %r12, 16(ARGa)
           
a3:         
   sub     ARGp, %r13
   jl      a4
   mov     %r13, 24(ARGa)
   
a4:
   pop     %r15
   pop     %r14
   pop     %r13
   pop     %r12
   pop     %r11
   pop     %r10
   pop     %rbx
   pop     %rax
   
#ifdef _WIN64 
   add     $120, %rsp
#endif

   ret