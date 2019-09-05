/* fpu_push_pop.S -- (C) Mark Rodenkirch, January 2018

   void      fpu_push_1divp(uint64_t p)
      Puts 1/p onto the FPU stack in st(0)
      
      Assumes FPU is set to double extended precision and round to zero.

   void      fpu_push_adivb(uint64_t a, uint64_t b)
      Puts a/b onto the FPU stack in st(0)
      
      Assumes FPU is set to double extended precision and round to zero.
      
   void      fpu_pop()
      Pop st(0)
     
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifdef _WIN64
#define ARGp      %rcx
#define ARGa      %rcx
#define ARGb      %rdx
#define tmp    32(%rsp)
#else
#define ARGp      %rdi
#define ARGa      %rdi
#define ARGb      %rsi
#define tmp   -32(%rsp)
#endif

   .p2align 4,,7
  
   .globl  _fpu_push_1divp
   .globl  fpu_push_1divp

_fpu_push_1divp:
fpu_push_1divp:
   mov     ARGp, tmp
  
   fildll  tmp
   fld1                      /* 1.0, p */
   fdivp   %st(0), %st(1)    /* 1.0/p, using round to zero. */

   ret

   
   
   .p2align 4,,7
  
   .globl  _fpu_push_adivb
   .globl  fpu_push_adivb

_fpu_push_adivb:
fpu_push_adivb:
   mov     ARGb, tmp
   fildll  tmp
   
   mov     ARGa, tmp
   fildll  tmp

   fdivp   %st(0), %st(1)    /* a/b, using round to zero. */

   ret


   
   .p2align 4,,7
  
  .globl  _fpu_pop
  .globl  fpu_pop

_fpu_pop:
fpu_pop:

   fstp    %st(0)

   ret
