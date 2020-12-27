/* inline.h -- (C) Mark Rodenkirch, June 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This is a collection of inline functions used by the framework.
*/

#ifndef _inline_H
#define _inline_H

inline uint32_t   gcd32(uint32_t a, uint32_t b)
{
   uint32_t c;

   while (b > 0)
   {
      c = a % b;
      a = b;
      b = c;
   }

   return a;
}

inline uint64_t  gcd64(uint64_t a, uint64_t b)
{
   uint64_t c;

   while (b > 0)
   {
      c = a % b;
      a = b;
      b = c;
   }

   return a;
}

// From gmp-fermat
inline uint64_t invproth(uint64_t k, uint32_t n)
{
   uint64_t i, j, x, y;

   x = (k<<n);
   y = ((k<<n)|1)<<n;

   for (i = 1, j = i<<n; j > 0; j <<= 1, y <<= 1)
      if (x & j)
         i += j, x += y;

   return -i;
}

// From gmp-fermat
inline uint64_t inv_mont(uint64_t y)
{
   uint64_t i, j, x;

   for (i = 1, j = 2, x = y-1; j > 0; j <<= 1)
   {
      y <<= 1;
      if (x & j)
         i += j, x += y;
   }

   return -i;
}

inline int32_t jacobi(int32_t a, uint32_t p)
{
   uint32_t x, y, t;
   int sign;

   if (a < 0)
      x = -a, sign = (p % 4 == 1) ? 1 : -1;
   else
      x = a, sign = 1;

   for (y = p; x > 0; x %= y)
   {
      for ( ; x % 2 == 0; x /= 2)
         if (y % 8 == 3 || y % 8 == 5)
            sign = -sign;

      t = x, x = y, y = t;

      if (x % 4 == 3 && y % 4 == 3)
         sign = -sign;
   }

   return ((y == 1) ? sign : 0);
}

// Thanks to Yves Gallot for this implementation based upon 
// Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
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

// Compute the residual of n (mod p)
inline uint64_t   mmmN(uint64_t n, uint64_t p)
{
   if (n == 1)
      return mmmOne(p);
    
   // list[0] = res(2^0), list[1] = res(2^1), list[2] = res(2^2), etc.
   uint64_t list[64];
   uint64_t bit = 0x01;
   uint64_t value = 0;

   list[0] = mmmOne(p);
   
   if (n & bit)
   {
      value = list[0];
      n &= ~bit;
   }
   
   for (uint32_t idx=1; idx<64; idx++)
   {
      bit <<= 1;
      
      // Need to compute for each power of 2
      list[idx] = mmmAdd(list[idx-1], list[idx-1], p);
      
      if (n & bit)
      {
         value = mmmAdd(list[idx], value, p);
         n &= ~bit;
      }
      
      if (n == 0)
         break;
   }

   return value;
}

// Compute the residual of b ^ n (mod p)
inline uint64_t   mmmPowmod(uint64_t base, uint64_t exp, uint64_t p, uint64_t q)
{
   uint64_t x = mmmN(base, p);
   uint64_t y = mmmN(1, p);
   
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
#endif
