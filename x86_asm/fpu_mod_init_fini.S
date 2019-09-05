/* fpu_mod_init_fini.S -- (C) Mark Rodenkirch, January 2018

   uint16_t  fpu_mod_init(void)
      Setup for multiplying modulo p, where 1 < p < 2^62.
      Sets FPU mode to double extended precision and round-to-zero.
      Returns old FPU mode precision and rounding bits.

   void      fpu_mod_fini(uint16_t mode)
      Restores old FPU mode precision and rounding bits.
     
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#define FPU_MODE_BITS (0x0F00)  /* 64-bit precision, round to zero. */

#ifdef _WIN64
#define ARG1   %cx
#define CW     8(%rsp)
#else
#define ARG1   %di
#define CW     -8(%rsp)
#endif

  .p2align 4,,7
  
  .globl  _fpu_mod_init
  .globl  fpu_mod_init

_fpu_mod_init:
fpu_mod_init:

   fnstcw  CW
   movzwl  CW, %eax            // Return old mode
   orw     $FPU_MODE_BITS, CW  // Set new mode bits
   fldcw   CW
   
   ret

  .p2align 4,,7
  
  .globl  _fpu_mod_fini
  .globl  fpu_mod_fini

_fpu_mod_fini:
fpu_mod_fini:

   fnstcw  CW
   and     $FPU_MODE_BITS, ARG1
   andw    $~FPU_MODE_BITS, CW // Clear bits
   or      ARG1, CW            // Restore old bits
   fldcw   CW

   ret
