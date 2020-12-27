/* mf_kernel.cl -- (C) Mark Rodenkirch, September 2016

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
 */

#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

void collect_factor(int  n,
                    int  c,
                    long p,
                    volatile __global int *factorsCount,
                     __global long4 *factors);

ulong mmmInvert(ulong p);
ulong mmmOne(ulong _p);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p);

__kernel void mf_kernel(__global const  ulong  *primes,
                        __global        ulong2 *rems,
                        __global        int    *params,
               volatile __global        int    *factorsCount,
                        __global        long4  *factors)
{
   int    gid = get_global_id(0);
   int    steps;
   int    idx = 0;
   int    maxNFirstLoop = D_MIN_N - D_MULTIFACTORIAL;
   
   ulong  startN = params[0];
   ulong  currentN = params[1];
   ulong  n;
   ulong  thePrime = primes[gid];

   if (thePrime == 0) return;

   ulong _q = mmmInvert(thePrime);
   ulong pOne = mmmOne(thePrime);
   ulong mOne = mmmSub(0, pOne, thePrime);

   ulong mfrs = mmmN(D_MULTIFACTORIAL, thePrime);;
   ulong ri, rf;
   
   if (startN == currentN)
   {
      // Set ri = residual of startN (mod p)
      // Set rf = residual of startN!mf (mod p)
      if (startN > thePrime)
         ri = mmmN(startN%thePrime, thePrime);
      else
         ri = mmmN(startN, thePrime);
      rf = ri;
   }
   else
   {      
      // Set ri = residual of currentN (mod p)
      // Set rf = residual of currentN!mf (mod p)
      ri = rems[gid].x;
      rf = rems[gid].y;
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
         collect_factor(n, -1, thePrime, factorsCount, factors);
         
      if (rf == mOne)
         collect_factor(n, +1, thePrime, factorsCount, factors);
   }

   rems[gid].x = ri;
   rems[gid].y = rf;
}

ulong mmmInvert(ulong p)
{
   ulong p_inv = 1;
   ulong prev = 0;
   
   while (p_inv != prev)
   {
      prev = p_inv;
      p_inv *= (2 - p * p_inv);
   }
   
   return p_inv;
}

// Compute the residual of 1 (mod p)
ulong mmmOne(ulong _p)
{
   return ((-_p) % _p);
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
   ulong lo = a * b;
   ulong hi = mul_hi(a, b);
   
   ulong m = lo * _q;
   
   ulong hi2 = mul_hi(m, _p);
   long r = (long) hi - (long) hi2;

   if (r < 0)
      return (ulong) (r + _p);
      
   return (ulong) r;
}

// Compute the residual of n (mod p)
ulong   mmmN(ulong n, ulong _p)
{
   if (n == 1)
      return mmmOne(_p);
    
   // list[0] = res(2^0), list[1] = res(2^1), list[2] = res(2^2), etc.
   ulong list[64];
   ulong bit = 0x01;
   ulong value = 0;

   list[0] = mmmOne(_p);
   
   if (n & bit)
   {
      value = list[0];
      n &= ~bit;
   }
   
   for (uint idx=1; idx<64; idx++)
   {
      bit <<= 1;
      
      // Need to compute for each power of 2
      list[idx] = mmmAdd(list[idx-1], list[idx-1], _p);
      
      if (n & bit)
      {
         value = mmmAdd(list[idx], value, _p);
         n &= ~bit;
      }
      
      if (n == 0)
         break;
   }

   return value;
}

void collect_factor(int  n,
                    int  c,
                    long p,
                    volatile __global int *factorsCount,
                     __global long4 *factors)
{
   if (n < D_MIN_N)
      return;
      
   if (n > D_MAX_N)
      return;
      
   int old = atomic_inc(factorsCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= D_MAX_FACTORS)
      return;
   
   factors[old].x = n;
   factors[old].y = c;
   factors[old].z = p;
}