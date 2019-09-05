/* asm-x86-fpu.h -- (C) Mark Rodenkirch, January 2018.

   Extended FPU x86 assembler routines for GCC.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   The 62-bit variant of the routines require p < 2^62.  To use the
   mulmod or powmod, you must call the mod_init to set the rounding
   mode and reset it with mod_fini when done.
*/

extern "C" {
   void      mod_52bit_fini(uint16_t mode);
   
   // Switch the FPU mode to double extended precision and round-to-zero
   uint16_t  fpu_mod_init(void);
   
   // Reset the FPU mode
   void      fpu_mod_fini(uint16_t mode);
   
   // Pushes 1/p onto top of FPU stack
   void      fpu_push_1divp(uint64_t p);
   
   // Pushes a/b onto top of FPU stack
   void      fpu_push_adivb(uint64_t a, uint64_t b);
   
   // Pop top of FPU stack
   void      fpu_pop(void);
   
   // Compute a*b (mod p).  Assumes that 1/p is at top of  FPU stack
   uint64_t  fpu_mulmod(uint64_t a, uint64_t b, uint64_t p);
   
   // Compute a*b (mod p).  Assumes b/p is at top of FPU stack.  Call this method
   // if you need compute each value between a*b^M (mod p) and a*b^N (mod p).
   uint64_t  fpu_mulmod_iter(uint64_t a, uint64_t b, uint64_t p);
   
   // Compute a*b (mod p) for 4 distinct a.  Assumes b/p is at top of FPU stack.
   // Call this method if you need compute each value between a*b^M (mod p) and
   // a*b^N (mod p) for 4 distinct a.
   void      fpu_mulmod_iter_4a(uint64_t *a, uint64_t b, uint64_t p);
   
   // Compute a*b (mod p) for 4 distinct a, b, and p.
   // Assumes FPU stack contains 1/p1, 1/p2, 1/p3, and 1/p4.
   void      fpu_mulmod_4a_4b_4p(uint64_t *a, uint64_t *b, uint64_t *p);
   
   // Compute b^n (mod p).  Assumes that 1/p is at top of FPU stack
   uint64_t  fpu_powmod(uint64_t b, uint64_t n, uint64_t p);
   
   // Compute b^n (mod p) for 4 distinct b and p.
   void      fpu_powmod_4b_1n_4p(uint64_t *b, uint64_t n, uint64_t *p);
}
