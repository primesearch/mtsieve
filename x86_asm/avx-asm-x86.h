/* asm-x86-avx.h -- (C) Mark Rodenkirch, January 2018.

   AVX x86 assembler routines for GCC.
   
   These routines require a, b, p < 2^52 and are designed to do
   powmod/mulmod/etc in groups of 16.
   
   They also require all double * parameters to be aligned to a 32-byte
   boundary such as shown here:
      double __attribute__((aligned(32))) p[16];

   Warning!  If the calling function or a function above is using a
             variable of type float or double, the contents of that variable
             could be changed by these routines as they do not save and
             restore the contents of the xmm/ymm/zmm registers.  Using
             double * is okay since those values are not stored in registers.
      
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#define AVX_ARRAY_SIZE  16

extern "C" {  
   // Given a list of 16 p, compute 1/p for each of them.
   void      avx_compute_reciprocal(double *p, double *reciprocal);
   
   // Set ymm12-ymm15 to the same value for a.
   void      avx_set_1a(double *a);
   
   // Set ymm12-ymm15 to the 16 different values for a.
   void      avx_set_16a(double *a);
   
   // Set ymm4-ymm7 to the 1 value for b.
   void      avx_set_1b(double *b);
   
   // Set ymm4-ymm7 to the 16 different values for b.
   void      avx_set_16b(double *b);
   
   // Set ymm12-ymm15 and b to b^n mod p.
   // b and p must be arrays of 16 doubles.
   void      avx_powmod(double *b, uint64_t n, double *p, double *reciprocal);
   
   // Set ymm12-ymm15 to a*b mod p.
   // Assumes that ymm12-ymm15 contain a.
   // Assumes that ymm4-ymm7 contain b.
   void      avx_mulmod(double *p, double *reciprocal);
   
   // Extract 64 double values from ymm0-ymm15
   void      avx_get_all(double *ptr);
   
   // Extract 16 double values from ymm12-ymm15
   void      avx_get_16a(double *ptr);
   
   // For the 16 doubles in ymm12-ymm15, compare to the comparator.
   // If a returned bit is 1, then the corresponding double in the 
   // ymm register is equal to the comparator.
   uint16_t  avx_pos_compare_1v(double *comparator);
   
   // Same as the 1v version, but 16 distinct comparators
   uint16_t  avx_pos_compare_16v(double *comparator);
   
   // For the 16 doubles in ymm12-ymm15, compare to the p - comparator.
   // If a returned bit is 1, then the corresponding double in the 
   // ymm register is equal to the comparator.
   uint16_t  avx_neg_compare_1v(double *comparator, double *p);
   
   // Same as the 1v version, but 16 distinct comparators
   uint16_t  avx_neg_compare_16v(double *comparator, double *p);
}
