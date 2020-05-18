/* fpu_powmod.S -- (C) Mark Rodenkirch, January 2018

   uint64_t fpu_powmod(uint64_t b, uint64_t n, uint64_t p);
      Returns b^n (mod p), where 0 <= b < p < 2^62.

      Assumes FPU is set to double extended precision and round to zero.
      Assumes that st(0) is 1/p.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifdef _WIN64
#define ARG1 %rcx
#define ARG2 %rdx
#define ARG3 %r8
#define REG4 %r9
#define REG5 %r10
#define REG6 %r11
#define TMP  16(%rsp) /* Shadow space */
#else
#define ARG1 %rdi
#define ARG2 %rsi
#define ARG3 %rdx
#define REG4 %rcx
#define REG5 %r8
#define REG6 %r9
#define TMP  -32(%rsp) /* Red zone */
#endif

   .text
   .globl _fpu_powmod
   .globl fpu_powmod

   .p2align 4,,15

_fpu_powmod:
fpu_powmod:

   push    %rbp
   
   lea     TMP, %rbp
   mov     $1, %eax
   mov     $1, REG4

mulsqr_loop:
   // %rax = REG4 = a, ARG1 = b, ARG2 = n, ARG3 = p, %st(0) = 1.0/p

   mov     ARG1, 0(%rbp)
   mov     %rax, 8(%rbp)
   
   fildll  0(%rbp)
   fildll  8(%rbp)
   
   fmul    %st(2), %st(0)
   fmul    %st(1), %st(0)
   fistpll 16(%rbp)
   
   fmul    %st(0), %st(0)
   fmul    %st(1), %st(0)
   fistpll 24(%rbp)
   
   imul    ARG1, %rax
   imul    ARG1, ARG1
   
   mov     16(%rbp), REG5
   mov     24(%rbp), REG6
   
   imul    ARG3, REG5
   imul    ARG3, REG6
   
   sub     REG5, %rax
   sub     REG6, ARG1
   
   mov     %rax, REG5
   mov     ARG1, REG6

correct_mul:
   sub     ARG3, REG5  /* CF=0 predicted */
   jl      correct_sqr
   mov     REG5, %rax

correct_sqr:
   sub     ARG3, REG6  /* CF=0 predicted */
   jl      shift
   mov     REG6, ARG1

shift:
   // REG4 = a, %rax = a*b, ARG1 = b^2, ARG2 = n

   shr     ARG2        /* CF unpredictable */
   cmovnc  REG4, %rax  /* Discard multiply if CF=0 */
   mov     %rax, REG4
   jnz     mulsqr_loop

   pop     %rbp
   
   ret
