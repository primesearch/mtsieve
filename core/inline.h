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

// Return (a/p) if gcd(a,p)==1, 0 otherwise.
inline int32_t jacobi(int64_t a, uint64_t p)
{
   uint64_t x, y, t;
   int32_t sign;

   if (a < 0)
   {
      x = -a;
      sign = (p % 4 == 1) ? 1 : -1;
   }
   else 
   {
      x = a;
      sign = 1;
   }

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

// Return the value of the Legendre symbol (a/p), where gcd(a,p)=1, p prime.
// If gcd(a,p)!=1 then return value undefined.
inline int32_t legendre(int64_t a, uint64_t p)
{
   uint64_t x, y, t;
   int32_t sign;

   if (a < 0)
   {
      a = -a;
      sign = (p % 4 == 1) ? 1 : -1;
   }
   else
      sign = 1;

   for (y = a; y % 2 == 0; y /= 2)
      if (p % 8 == 3 || p % 8 == 5)
         sign = -sign;

   if (p % 4 == 3 && y % 4 == 3)
      sign = -sign;

   for (x = p % y; x > 0; x %= y)
   {
      for ( ; x % 2 == 0; x /= 2)
         if (y % 8 == 3 || y % 8 == 5)
            sign = -sign;

      t = x, x = y, y = t;

      if (x % 4 == 3 && y % 4 == 3)
         sign = -sign;
   }

   return sign;
}

#endif
