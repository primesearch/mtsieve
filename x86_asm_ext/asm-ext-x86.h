/* asm-ext-x86.h -- (C) Mark Rodenkirch, June 2018.

   Extended x86 assembler routines for GCC.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   The 62-bit variant of the routines require p < 2^62.  To use the
   mulmod or powmod, you must call the mod_init to set the rounding
   mode and reset it with mod_fini when done.
*/

// Handle the possibility that mp_limb_t and uint64_t are different typedefs.
// We don't use mp_limb_t here, but the same parameters passed to these calls
// are also passed to GMP functions so the types must be compatible.
#ifdef WIN32
   #define uns_64 unsigned long long
#else
   #define uns_64 unsigned long
#endif

extern "C" {
   // Assembly language functions in mulmod128.S
   extern void mulmod128(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N, const uns_64 inv);
   extern void mulmod128_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in sqrmod128.S
   extern void sqrmod128(const uns_64 *R, const uns_64 *N, const uns_64 inv);
   extern void sqrmod128_proth0(const uns_64 *R, const uns_64 *N, const uns_64 inv);
   extern void sqrmod128_proth1(const uns_64 *R, const uns_64 *N);

   // Assembly language functions in mulmod192.S
   extern void mulmod192(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N, const uns_64 inv);
   extern void mulmod192_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod192_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in sqrmod192.S
   extern void sqrmod192(const uns_64 *R, const uns_64 *N, const uns_64 inv);
   extern void sqrmod192_proth0(const uns_64 *R, const uns_64 *N);
   extern void sqrmod192_proth1(const uns_64 *R, const uns_64 *N);

   // Assembly language functions in mulmod256.S
   extern void mulmod256_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod256_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in sqrmod256.S
   extern void sqrmod256_proth0(const uns_64 *R, const uns_64 *N);
   extern void sqrmod256_proth1(const uns_64 *R, const uns_64 *N);

   // Assembly language functions in m320.S
   extern void mulmod320_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod320_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m384.S
   extern void mulmod384_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod384_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m448.S
   extern void mulmod448_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod448_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m512.S
   extern void mulmod512_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod512_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m576.S
   extern void mulmod576_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod576_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m640.S
   extern void mulmod640_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod640_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m704.S
   extern void mulmod704_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod704_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in m768.S
   extern void mulmod768_proth0(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);
   extern void mulmod768_proth1(const uns_64 *R, const uns_64 *A, const uns_64 *B, const uns_64 *N);

   // Assembly language functions in redc.S
   extern void redc_proth0(uns_64 *R, const uns_64 *N, size_t nn);
}
