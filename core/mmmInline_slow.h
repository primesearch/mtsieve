/* mmmInline.h -- (C) Mark Rodenkirch, June 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This is a collection of inline functions to are used for
   Montgomery multiplication used by the framework.

   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
*/

#ifndef _inline_H
#define _inline_H

inline uint64_t   mmmInvert(uint64_t p)
{
   return 0;
}

// Compute the residual of 1 (mod p)
inline uint64_t   mmmOne(uint64_t p)
{
   return 1;
}

// Compute the residual of a + b (mod p)
// Note that resA is residual of a (mod p) and resB is residual of b (mod p)
inline uint64_t   mmmAdd(uint64_t resA, uint64_t resB, uint64_t p)
{
   uint64_t c = (resA >= p - resB) ? p : 0;
   return resA + resB - c;
}

// Compute the residual of a - b (mod p)
// Note that resA is residual of a (mod p) and resB is residual of b (mod p)
inline uint64_t   mmmSub(uint64_t resA, uint64_t resB, uint64_t p)
{
   uint64_t c = (resA < resB) ? p : 0;
   return resA - resB + c;
}

// Compute the residual of a * b (mod p)
// Note that resA is residual of a (mod p) and resB is residual of b (mod p)
inline uint64_t   mmmMulmod(uint64_t resA, uint64_t resB, uint64_t p, uint64_t q)
{
   return uint64_t((resA * __uint128_t(resB)) % p);
}

// Compute the residual of b ^ n (mod p)
// Note that resB is residual of b (mod p) and resOne is residual of 1 (mod p)
inline uint64_t   mmmPowmod(uint64_t resB, uint64_t exp, uint64_t p, uint64_t q, uint64_t resOne)
{
   uint64_t x = resB;
   uint64_t y = resOne;
   
   while (true)
   {
      if (exp & 1)
         y = mmmMulmod(x, y, p, q);

      exp >>= 1;

      if (!exp)
         return y;

      x = mmmMulmod(x, x, p, q);
   }

   // Should never get here
   return 0;
}

// Set b to the residual of b ^ n (mod p)
// Note that resB is residual of b (mod p) and resOne is residual of 1 (mod p)
// resBexpP is output only and the other parameters are input only.
inline void   mmmPowmodX4(uint64_t *resB, uint64_t exp, uint64_t *p, uint64_t *q, uint64_t *resOne, uint64_t *resBexpP)
{
   uint64_t resX0 = resB[0];
   uint64_t resX1 = resB[1];
   uint64_t resX2 = resB[2];
   uint64_t resX3 = resB[3];
   
   uint64_t resY0 = resOne[0];
   uint64_t resY1 = resOne[1];
   uint64_t resY2 = resOne[2];
   uint64_t resY3 = resOne[3];
   
   while (true)
   {
      if (exp & 1)
      {
         resY0 = mmmMulmod(resX0, resY0, p[0], q[0]);
         resY1 = mmmMulmod(resX1, resY1, p[1], q[1]);
         resY2 = mmmMulmod(resX2, resY2, p[2], q[2]);
         resY3 = mmmMulmod(resX3, resY3, p[3], q[3]);
      }

      exp >>= 1;

      if (!exp)
         break;

      resX0 = mmmMulmod(resX0, resX0, p[0], q[0]);
      resX1 = mmmMulmod(resX1, resX1, p[1], q[1]);
      resX2 = mmmMulmod(resX2, resX2, p[2], q[2]);
      resX3 = mmmMulmod(resX3, resX3, p[3], q[3]);
   }
   
   resBexpP[0] = resY0;
   resBexpP[1] = resY1;
   resBexpP[2] = resY2;
   resBexpP[3] = resY3;
}

// Compute the residual of 2^64 (mod p)
inline uint64_t   mmm2exp64(uint64_t p, uint64_t q, uint64_t resOne)
{
   return 0;
}

// Compute the residual of n (mod p)
inline uint64_t   mmmN(uint64_t n, uint64_t p, uint64_t q, uint64_t res2exp64)
{
   return n;
}

#endif
