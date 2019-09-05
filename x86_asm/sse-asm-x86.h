/* asm-x86-sse.h -- (C) Mark Rodenkirch, January 2018.

   SSE x86 assembler routines for GCC.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   The 52-bit variant of the routines require p < 2^52.  To use the
   mulmod or powmod, you must call the mod_init to set the rounding
   mode and reset it with mod_fini when done.
*/

extern "C" {

   // Switch the SSE mode to round-to-zero
   uint16_t  sse_mod_init(void);
   
   // Reset the SSE mode
   void      sse_mod_fini(uint16_t mode);
   
   // Compute a = a*b (mod p) for given a, b, and p.
   uint64_t  sse_mulmod(uint64_t a, uint64_t b, uint64_t p, double *bdivp);
   
   // Compute a = a*b (mod p) for 4 distinct a, b, and p.
   void      sse_mulmod_4a_4b_4p(uint64_t *a, uint64_t *b, uint64_t *p, double *invp);
   
   // Compute b = b^n (mod p) for 4 distinct b and p.
   void      sse_powmod_4b_1n_4p(uint64_t *b, uint64_t n, uint64_t *p, double *invp);
   
   // Compute b = mul*b^n (mod p) for 4 distinct b and p.
   void      sse_powmod_4b_1n_4p_mulmod_1k(uint64_t *b, uint64_t n, uint64_t *p, uint64_t mul);
}
