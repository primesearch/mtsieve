/* mf_kernel.metal -- (C) Mark Rodenkirch, April 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
 */

#include <metal_stdlib>
#include <metal_atomic>
using namespace metal;

void collect_factor(int         n,
                    int         c,
                    long        p,
             device int        *params,
    volatile device atomic_int *factorsCount,
             device long4      *factors);

ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);

#define P_MIN_N            0
#define P_MAX_N            1
#define P_MAX_STEPS        2
#define P_MULTIFACTORIAL   3
#define P_MAX_FACTORS      4
#define P_START_N          5
#define P_CURRENT_N        6

kernel void mf_kernel(device const  ulong  *primes,
                      device        ulong2 *rems,
                      device        int    *params,
             volatile device atomic_int    *factorsCount,
                      device        long4  *factors,
                                    uint    tpid [[thread_position_in_grid]])
{
   int    steps;
   ulong  D_MIN_N = params[P_MIN_N];
   ulong  D_MAX_N = params[P_MAX_N];
   int    D_MULTIFACTORIAL = params[P_MULTIFACTORIAL];
   int    D_MAX_STEPS = params[P_MAX_STEPS];
   ulong  maxNFirstLoop = D_MIN_N - D_MULTIFACTORIAL;
   
   ulong  startN = params[P_START_N];
   ulong  currentN = params[P_CURRENT_N];
   ulong  n;
   ulong  thePrime = primes[tpid];

   if (thePrime == 0) return;

   ulong _q = mmmInvert(thePrime);
   ulong pOne = mmmOne(thePrime);
   ulong mOne = mmmSub(0, pOne, thePrime);
   ulong _r2 = mmmR2(thePrime, _q, pOne);

   ulong mfrs = mmmN(D_MULTIFACTORIAL, thePrime, _q, _r2);
   ulong ri, rf;
   
   if (startN == currentN)
   {
      // Set ri = residual of startN (mod p)
      // Set rf = residual of startN!mf (mod p)
      if (startN > thePrime)
         ri = mmmN(startN%thePrime, thePrime, _q, _r2);
      else
         ri = mmmN(startN, thePrime, _q, _r2);
      rf = ri;
   }
   else
   {      
      // Set ri = residual of currentN (mod p)
      // Set rf = residual of currentN!mf (mod p)
      ri = rems[tpid].x;
      rf = rems[tpid].y;
   }

   steps = 0;
   
   // ri = residual of currentN (mod p)
   // rf = residual of currentN!mf (mod p)
   
   n = currentN + D_MULTIFACTORIAL;
   
   for (; n<maxNFirstLoop; n+=D_MULTIFACTORIAL, steps++)
   {
      if (steps >= D_MAX_STEPS)
         break;
         
      ri = mmmAdd(ri, mfrs, thePrime);
      rf = mmmMulmod(rf, ri, thePrime, _q);
   }

   for (; n <=D_MAX_N; n+=D_MULTIFACTORIAL, steps++)   
   {      
      if (steps >= D_MAX_STEPS)
         break;

      ri = mmmAdd(ri, mfrs, thePrime);
      rf = mmmMulmod(rf, ri, thePrime, _q);

      if (rf == pOne)
         collect_factor(n, -1, thePrime, params, factorsCount, factors);
         
      if (rf == mOne)
         collect_factor(n, +1, thePrime, params, factorsCount, factors);
   }

   rems[tpid].x = ri;
   rems[tpid].y = rf;
}

ulong mmmInvert(ulong _p)
{
   ulong p_inv = 1;
   ulong prev = 0;
   
   while (p_inv != prev)
   {
      prev = p_inv;
      p_inv *= (2 - _p * p_inv);
   }
   
   return p_inv;
}

// Compute the residual of 1 (mod p)
ulong mmmOne(ulong _p)
{
   return ((-_p) % _p);
}

// Compute the residual of 2^64 (mod p)
ulong mmmR2(ulong _p, ulong _q, ulong _one)
{
	ulong t = mmmAdd(_one, _one, _p);
   
   t = mmmAdd(t, t, _p);   // 4
	for (size_t i=0; i<5; i++)
      t = mmmMulmod(t, t, _p, _q);   // 4^{2^5} = 2^64
      
	return t;
}

ulong mmmAdd(ulong a, ulong b, ulong _p)
{
   ulong c = (a >= _p - b) ? _p : 0;
   return a + b - c;
}

ulong mmmSub(ulong a, ulong b, ulong _p)
{
   ulong c = (a < b) ? _p : 0;
   return a - b + c;
}

ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q)
{
   ulong lo, hi;

   lo = a * b; hi = mulhi(a, b);

   ulong m = lo * _q;
   ulong mp = mulhi(m, _p);
   long r = (long)(hi - mp);
   return (r < 0) ? r + _p : r;
}

// Compute the residual of n (mod p)
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2)
{
   return mmmMulmod(n, _r2, _p, _q);
}

void collect_factor(int         n,
                    int         c,
                    long        p,
             device int        *params,
    volatile device atomic_int *factorsCount,
             device long4      *factors)
{
   int    D_MIN_N = params[P_MIN_N];
   int    D_MAX_N = params[P_MAX_N];
   int    D_MAX_FACTORS = params[P_MAX_FACTORS];
   
   if (n < D_MIN_N)
      return;
      
   if (n > D_MAX_N)
      return;
      
   int old = atomic_fetch_add_explicit(factorsCount, 1, memory_order_relaxed);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= D_MAX_FACTORS)
      return;
   
   factors[old].x = n;
   factors[old].y = c;
   factors[old].z = p;
}
