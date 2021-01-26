/* MpArith.h -- (C) Mark Rodenkirch, January 2021

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
*/

#ifndef _MpArith_H
#define _MpArith_H

//#define DEBUG_MP

typedef uint64_t  MpRes;

class MpArith
{
private:
   const uint64_t _p;
   const uint64_t _q;
   const MpRes _one; 
   MpRes _r2exp64;

private:
   // p * p_inv = 1 (mod 2^64) (Newton's method)
   static uint64_t invert(const uint64_t p)
   {
#ifdef DEBUG_MP
      return 0;
#else
      uint64_t p_inv = 1, prev = 0;

      while (p_inv != prev)
      {
         prev = p_inv;
         p_inv *= 2 - p * p_inv;
      }

      return p_inv;
#endif
   }

   uint64_t REDC(const __uint128_t t) const
   {
#ifdef DEBUG_MP
      return uint64_t(t % _p);
#else
      const uint64_t m = uint64_t(t) * _q;
      const int64_t r = int64_t((t >> 64) - uint64_t((m * __uint128_t(_p)) >> 64));

      return (r < 0) ? uint64_t(r + _p) : uint64_t(r);
#endif
   }

public:

#ifdef DEBUG_MP
   MpArith(const uint64_t p) : _p(p), _q(invert(p)), _one(1)
#else
   MpArith(const uint64_t p) : _p(p), _q(invert(p)), _one((-p) % p)
#endif
   {
      MpRes t = add(_one, _one);
      
      t = add(t, t); // 4
		for (size_t i = 0; i < 5; ++i)
         t = mul(t, t); // 4^{2^5} = 2^64
         
		_r2exp64 = t;
   } 

   uint64_t p() const { return _p; }
   
   MpRes one() const { return _one; }

   MpRes add(const MpRes a, const MpRes b) const
   {
      MpRes r;
      
      const uint64_t c = (a >= _p - b) ? _p : 0;

      r = a + b - c;

      return r;
   }

   MpRes sub(const MpRes a, const MpRes b) const
   {
      MpRes r;
      
      const uint64_t c = (a < b) ? _p : 0;

      r = a - b + c;

      return r;
   }

   MpRes mul(const MpRes a, const MpRes b) const
   {
      MpRes r;
      
      r = REDC(a * __uint128_t(b));

      return r;
   }

	MpRes pow(MpRes a, size_t exp) const
	{
      uint64_t x = a;
      uint64_t y = _one;
      MpRes r;

      while (true)
      {
         if (exp & 1)
            y = REDC(x * __uint128_t(y));

         exp >>= 1;

         if (!exp)
            break;

         x = REDC(x * __uint128_t(x));
      }
      
      r = y;
      return r;
	}
   
   // Compute the residual of n (mod p)
   MpRes nToRes(uint64_t n)
   {
#ifdef DEBUG_MP
      return n;
#else
      if (n == 1)
         return _one;
       
      return mul(n, _r2exp64);
#endif
   }

   uint64_t resToN(MpRes res)
   {
#ifdef DEBUG_MP
      return res;
#else
      MpRes r;
      
      r = REDC(res);
      
      return r;
#endif
   }
};

#endif