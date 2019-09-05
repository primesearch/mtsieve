/* sse_mod_init_fini.S -- (C) Mark Rodenkirch, May 2018

   uint16_t  sse_mod_init(void)
      Setup for multiplying modulo p, where 1 < p < 2^52.
      Sets SSE mode to round-to-zero.
      Returns old SSE mode rounding bits.

   void      sse_mod_fini(uint16_t mode)
      Restores old SSE mode rounding bits.
     
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#define SSE_MODE_BITS (0x6000)		/* Round To Zero */

#ifdef _WIN64
#define ARG1   %ecx

#define OLD_MXCSR  4(%rsp)
#define NEW_MXCSR  8(%rsp)
#else
#define ARG1   %edi

#define OLD_MXCSR -4(%rsp)
#define NEW_MXCSR -8(%rsp)
#endif

  .p2align 4,,7
  
  .globl  _sse_mod_init
  .globl  sse_mod_init

_sse_mod_init:
sse_mod_init:

#ifdef _WIN64
   sub      $16, %rsp
#endif

   stmxcsr OLD_MXCSR
   mov     OLD_MXCSR, %eax
   or      $SSE_MODE_BITS, %eax      /* Round to zero */
   mov     %eax, NEW_MXCSR
   ldmxcsr NEW_MXCSR
   mov     OLD_MXCSR, %eax
   
#ifdef _WIN64
   add      $16, %rsp
#endif

   ret

  .p2align 4,,7
  
  .globl  _sse_mod_fini
  .globl  sse_mod_fini

_sse_mod_fini:
sse_mod_fini:

#ifdef _WIN64
   sub      $16, %rsp
#endif

   stmxcsr OLD_MXCSR
   mov	  OLD_MXCSR, %eax
   mov     ARG1, NEW_MXCSR
   and	  $~SSE_MODE_BITS, %eax
   or      %eax, NEW_MXCSR 	      /* Restore old bits */
   ldmxcsr NEW_MXCSR

#ifdef _WIN64
   add      $16, %rsp
#endif

   ret
