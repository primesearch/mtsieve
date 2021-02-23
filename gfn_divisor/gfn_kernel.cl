/* gfn_kernel.cl -- (C) Mark Rodenkirch, December 2013

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

void  collectFactor(ulong k, uint n, ulong p, 
     volatile __global       uint   *factorCount,
              __global       ulong4 *factors);
              
uint  getSmallDivisor(ulong k, uint n);
ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);
ulong mmmResToN(ulong res, ulong _p, ulong _q);

// Note that primes[gid] > K_MAX is a requirement when using these kernels

#ifndef D_MULTI_PASS
__kernel void gfn_kernel(__global const ulong  *primes,
                volatile __global       uint   *factorCount,
                         __global       ulong4 *factors)
{
   int    gid = get_global_id(0);
   uint   n;
   ulong  p = primes[gid];
   ulong  k = (1+p) >> 1;
   ushort bitsToShift;

   ulong _q = mmmInvert(p);
   ulong _one = mmmOne(p);
   ulong _r2 = mmmR2(p, _q, _one);
   
   ulong mpK = mmmNtoRes(k, p, _q, _r2);

   ulong mpRem = mmmPowmod(mpK, N_MIN, p, _q, _one, _r2);

   k = p - mmmResToN(mpRem, p, _q);
   n = N_MIN;

   while (n <= N_MAX)
   {
      // Determine how many times we have to divide k by 2
      // to ensure that k is odd.
      bitsToShift = 63 - clz(k & -k);
      
      k >>= bitsToShift;
      n += bitsToShift;
      
      if (k <= K_MAX && k >= K_MIN && n <= N_MAX)
         collectFactor(k, n, p, factorCount, factors);

      // Make k even before the next iteration of the loop
      k += p;
   }
}

#else
// Note that primes[gid] > K_MAX is a requirement when using this kernel
__kernel void gfn_kernel(__global const ulong  *primes,
                         __global const ulong  *params,
                         __global       ulong  *rems,
                volatile __global       uint   *factorCount,
                         __global       ulong4 *factors)
{
   int    gid = get_global_id(0);
   uint   n, steps = 0;
   ulong  p = primes[gid];
   ulong  k = (1+p) >> 1;
   ushort bitsToShift;

   if (params[0] == N_MIN)
   {
      ulong _q = mmmInvert(p);
      ulong _one = mmmOne(p);
      ulong _r2 = mmmR2(p, _q, _one);
      
      ulong mpK = mmmNtoRes(k, p, _q, _r2);

      ulong mpRem = mmmPowmod(mpK, N_MIN, p, _q, _one, _r2);

      k = p - mmmResToN(mpRem, p, _q);
      n = N_MIN;
   }
   else
   {
      k = rems[gid];
      n = params[0];
   }

   while (n <= N_MAX - 64 && steps < D_MAX_STEPS - 64)
   {
      // Determine how many times we have to divide k by 2
      // to ensure that k is odd.
      bitsToShift = 63 - clz(k & -k);
      
      k >>= bitsToShift;
      n += bitsToShift;
      
      steps += bitsToShift;
      
      if (k <= K_MAX && k >= K_MIN)
         collectFactor(k, n, p, factorCount, factors);

      // Make k even before the next iteration of the loop
      k += p;
   }

   while (n <= N_MAX && steps < D_MAX_STEPS)
   {
      // We only need to collect odd k
      if (k & 1)
      {
         if (k <= K_MAX && k >= K_MIN && n <= N_MAX)
            collectFactor(k, n, p, factorCount, factors);

         // Make k even before dividing by 2
         k += p;
      }

      // For the next n, divide k by 2
      // Note that k*2^n+1 = (k/2)*2^(n+1)+1

      k >>= 1;
      n++;
      steps++;
   }

   rems[gid] = k;
}
#endif

void  collectFactor(ulong k, uint n, ulong p, 
     volatile __global       uint   *factorCount,
              __global       ulong4 *factors)
{
   // Do not report k/n if it has a known small divisor
   if (getSmallDivisor(k, n) > 0)
      return;

   int old = atomic_inc(factorCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= D_MAX_FACTORS)
      return;

   factors[old].x = k;
   factors[old].y = n;
   factors[old].z = p;
}

// Check for small divisors.
uint  getSmallDivisor(ulong k, uint n)
{
   uint  smallN;
   ulong smallK;
   
   smallN = n % (2);
   smallK = k % (3);
   if ((smallK << smallN) % (3) == 2) return 3;
   
   smallN = n % (4);
   smallK = k % (5);
   if ((smallK << smallN) % (5) == 4) return 5;
   
   smallN = n % (6);
   smallK = k % (7);
   if ((smallK << smallN) % (7) == 6) return 7;
   
   smallN = n % (10);
   smallK = k % (11);
   if ((smallK << smallN) % (11) == 10) return 11;
   
   smallN = n % (12);
   smallK = k % (13);
   if ((smallK << smallN) % (13) == 12) return 13;
   
   smallN = n % (16);
   smallK = k % (17);
   if ((smallK << smallN) % (17) == 16) return 17;
   
   smallN = n % (18);
   smallK = k % (19);
   if ((smallK << smallN) % (19) == 18) return 19;
   
   smallN = n % (22);
   smallK = k % (23);
   if ((smallK << smallN) % (23) == 22) return 23;
   
   smallN = n % (28);
   smallK = k % (29);
   if ((smallK << smallN) % (29) == 28) return 29;
   
   smallN = n % (30);
   smallK = k % (31);
   if ((smallK << smallN) % (31) == 30) return 31;

   smallN = n % (36);
   smallK = k % (37);
   if ((smallK << smallN) % (37) == 36) return 37;
   
   smallN = n % (40);
   smallK = k % (41);
   if ((smallK << smallN) % (41) == 40) return 41;
   
   smallN = n % (42);
   smallK = k % (43);
   if ((smallK << smallN) % (43) == 42) return 43;
   
   smallN = n % (46);
   smallK = k % (47);
   if ((smallK << smallN) % (47) == 46) return 47;
   
   return 0;
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
ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2)
{
   return mmmMulmod(n, _r2, _p, _q);
}

// Convert a residual back to n
ulong mmmResToN(ulong res, ulong _p, ulong _q)
{
   return mmmMulmod(res, 1, _p, _q);
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

