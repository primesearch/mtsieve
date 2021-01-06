/* gfn_kernel.cl -- (C) Mark Rodenkirch, December 2013

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

void collectFactor(ulong   k,
                   uint    n,
                   ulong   p,
 volatile __global uint   *factorCount,
          __global ulong4 *factors);
          
ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);

// Note that primes[gid] > K_MAX is a requirement when using this kernel
__kernel void gfn_kernel(__global const ulong  *primes,
                volatile __global       uint   *factorCount,
                         __global       ulong4 *factors)
{
   int    gid = get_global_id(0);
   uint   n;
   ulong  p = primes[gid];
   ulong  k = (1+p) >> 1;

   ulong _q = mmmInvert(p);
   ulong _one = mmmOne(p);
   ulong _r2 = mmmR2(p, _q, _one);
   
   ulong mpK = mmmN(k, p, _q, _r2);

   ulong mpRes = mmmPowmod(mpK, N_MIN, p, _q, _one, _r2);

   ulong rem = mmmMulmod(mpRes, 1, p, _q);
      
   k = p - rem;

   for (n=N_MIN; n<=N_MAX; n++)
   {
      if (k & 1)
      {
         // We only need to collect even k
         if (k <= K_MAX && k >= K_MIN)
           collectFactor(k, n, p, factorCount, factors);

         k += p;
      }
      
      // For the next n, divide k by 2
      // Note that k*2^n+1 = (k/2)*2^(n+1)+1

      k >>= 1;
   }
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

// Compute the residual of b ^ n (mod p)
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2)
{
   ulong x = resB;
   ulong y = _one;

   while (true)
   {
      if (exp & 1)
         y = mmmMulmod(x, y, _p, _q);

      exp >>= 1;

      if (!exp)
         return y;

      x = mmmMulmod(x, x, _p, _q);
   }

   // Should never get here
   return 0;
}

void collectFactor(ulong   k,
                   uint    n,
                   ulong   p,
 volatile __global uint   *factorCount,
          __global ulong4 *factors)
{
   int old = atomic_inc(factorCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= MAX_FACTORS)
      return;

   factors[old].x = k;
   factors[old].y = n;
   factors[old].z = p;
}

