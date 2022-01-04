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

ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);

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

#ifdef __NV_CL_C_VERSION
   const uint a0 = (uint)(a), a1 = (uint)(a >> 32);
   const uint b0 = (uint)(b), b1 = (uint)(b >> 32);

   uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;

   asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
   asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c2) : "r" (a0), "r" (b1));

   asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b1), "r" (c2));
   asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c3) : "r" (a1), "r" (b1));

   asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
   asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
   asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c3) : "r" (c3));

   lo = upsample(c1, c0); hi = upsample(c3, c2);
#else
   lo = a * b; hi = mul_hi(a, b);
#endif

   ulong m = lo * _q;
   ulong mp = mul_hi(m, _p);
   long r = (long)(hi - mp);
   return (r < 0) ? r + _p : r;
}

// Compute the residual of n (mod p)
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2)
{
   return mmmMulmod(n, _r2, _p, _q);
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