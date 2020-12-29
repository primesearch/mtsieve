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
   uint64_t p_inv = 1;
   uint64_t prev = 0;
   
   while (p_inv != prev)
   {
      prev = p_inv;
      p_inv *= (2 - p * p_inv);
   }
   
   return p_inv;
}

// Compute the residual of 1 (mod p)
inline uint64_t   mmmOne(uint64_t p)
{
   return ((-p) % p);
}

// Compute the residual of a + b (mod p)
inline uint64_t   mmmAdd(uint64_t a, uint64_t b, uint64_t p)
{
   uint64_t c = (a >= p - b) ? p : 0;
   return a + b - c;
}

// Compute the residual of a - b (mod p)
inline uint64_t   mmmSub(uint64_t a, uint64_t b, uint64_t p)
{
   uint64_t c = (a < b) ? p : 0;
   return a - b + c;
}

// Compute the residual of a * b (mod p)
inline uint64_t   mmmMulmod(uint64_t a, uint64_t b, uint64_t p, uint64_t q)
{
   __uint128_t t1 = a * __uint128_t(b);
   
   uint64_t m = uint64_t(t1) * q;
   
   __uint128_t t2 = m * __uint128_t(p);
   
   int64_t r = int64_t(t1 >> 64) - int64_t(t2 >> 64);

   if (r < 0)
      return (uint64_t) (r + p);
      
   return (uint64_t) r;      
}
// Compute the residual of b ^ n (mod p)
inline uint64_t   mmmPowmod(uint64_t resB, uint64_t exp, uint64_t p, uint64_t q)
{
   uint64_t x = resB;
   uint64_t y = mmmOne(p);
   
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

// Compute the residual of 2^64 (mod p)
inline uint64_t   mmm2exp64(uint64_t p, uint64_t q, uint64_t res1)
{
	uint64_t t = mmmAdd(res1, res1, p);
   
   // t = res 4 (mod p)
   t = mmmAdd(t, t, p);   // 4
   
	for (uint32_t i=0; i<5; i++)
      t = mmmMulmod(t, t, p, q);   // 4^{2^5} = 2^64
      
	// t = res 2^64 (mod p)
   return t;
}

// Compute the residual of n (mod p)
inline uint64_t   mmmN(uint64_t n, uint64_t p, uint64_t q, uint64_t res2exp64)
{
   if (n == 1)
      return mmmOne(p);
    
   return mmmMulmod(n, res2exp64, p, q);
}

#endif
