/* cw_kernel.cl -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

void collect_factor(int  n,
                    int  c,
                    long p,
                    volatile __global uint *factorCount,
                    __global long4 *factors);

long compute_inverse(long a, long p);

ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);

#define MAX_POWERS 50

__kernel void cw_kernel(__global const  long  *primes,
                        __global const  int   *terms,
               volatile __global        uint  *factorCount,
                        __global        long4 *factors)
{
   int gid = get_global_id(0);

   ulong thePrime = primes[gid];

   ulong powers[MAX_POWERS+1];
   ulong resC[MAX_POWERS+1];
   ulong resW[MAX_POWERS+1];

   // We're trying to determine if n*b^n (mod p) = +1/-1.  This requires an two
   // operations, an exponentiation (b^n) and a multiplcation (n).  We can make a
   // change to eliminate one of those operations.  Doing this requires us to
   // compute the multiplicative inverse of b.  The multiplicated inverse of
   // b (which I will call B) is the number such at b*B = 1 mod p.  Here is how
   // the algorithm works.
   //
   // We start with our number:
   //       n*b^n (mod p) = +1/-1
   //    -> n*b^n*B^n (mod p) = +1/-1 * B^n (mod p)
   //    -> n (mod p) = +1/-1 B^n (mod p)
   //    -> +n/-n (mod p) = B^n (mod p)
   // And voila, if B^n (mod p) = -n, then p is a factor of n*b^n+1.
   // and if B^n (mod p) = +n, then p is a factor of n*b^n-1.
   //
   // To get then next term (which I will call N) where n > N, we need to
   // multiply B^n by B^(N-n):
   //    -> +n/-n (mod p) = B^n * B(N-n) (mod p)
   // Since N-n is negative we could compute (p-1)+N-n, but that could be
   // a very large exponent.  Instead since B is the multiplicative inverse
   // of b, we can use that, i.e. b^(n-N), which will have a much smaller
   // exponent and that allows us to put b^(n-N) into a table so that we
   // can use mulmods instead of expmods as we iterate through n.
   //    -> +n/-n (mod p) = B^n * b^n (mod p)

   ulong _q = mmmInvert(thePrime);
   ulong _one = mmmOne(thePrime);
   ulong _r2 = mmmR2(thePrime, _q, _one);

   // compute the inverse of b (mod p)
   ulong inverse = compute_inverse(BASE, thePrime);

   int theN = terms[0];

   ulong resI = mmmN(inverse, thePrime, _q, _r2);
   ulong rem = mmmPowmod(resI, theN, thePrime, _q, _one, _r2);

#ifdef CHECK_CULLEN
   ulong cullen  = mmmN(thePrime - theN, thePrime, _q, _r2);

   if (rem == cullen)
      collect_factor(theN, +1, thePrime, factorCount, factors);
#endif

#ifdef CHECK_WOODALL
   ulong woodall = mmmN(theN, thePrime, _q, _r2);

   if (rem == woodall)
      collect_factor(theN, -1, thePrime, factorCount, factors);
#endif

   powers[1] = mmmN(BASE, thePrime, _q, _r2);
   powers[2] = mmmMulmod(powers[1], powers[1], thePrime, _q);

   resC[1] = mmmN(thePrime - 1, thePrime, _q, _r2);
   resW[1] = _one;

   resC[2] = mmmAdd(resC[1], resC[1], thePrime);
   resW[2] = mmmAdd(resW[1], resW[1], thePrime);

   if (BASE & 1)
   {
      // If the base is odd, then all n must be even and thus the difference
      // between any two remaining n for the base must also be even.
      for (size_t idx=4; idx<MAX_POWERS+1; idx+=2)
      {
         powers[idx] = mmmMulmod(powers[idx-2], powers[2], thePrime, _q);
         resC[idx] = mmmAdd(resC[idx-2], resC[2], thePrime);
         resW[idx] = mmmAdd(resW[idx-2], resW[2], thePrime);
      }
   }
   else
   {
      for (size_t idx=3; idx<MAX_POWERS+1; idx++)
      {
        powers[idx] = mmmMulmod(powers[idx-1], powers[1], thePrime, _q);
        resC[idx] = mmmAdd(resC[idx-1], resC[1], thePrime);
        resW[idx] = mmmAdd(resW[idx-1], resW[1], thePrime);
      }
   }

   int prevN = terms[0];

   // Note that the terms start at max N and decrease
   size_t termIndex = 1;
   while (terms[termIndex] > 0)
   {
      theN = terms[termIndex];

      size_t power = (size_t)(prevN - theN);

      // Our table isn't infinite in size, so we'll break down the power into smaller
      // chunks.  At worst this might cause an extra mulmod or two every once in a while.
      while (power > MAX_POWERS)
      {
         rem = mmmMulmod(rem, powers[MAX_POWERS], thePrime, _q);
#ifdef CHECK_CULLEN
         cullen = mmmSub(cullen, resC[MAX_POWERS], thePrime);
#endif

#ifdef CHECK_WOODALL
         woodall = mmmSub(woodall, resW[MAX_POWERS], thePrime);
#endif
         power -= MAX_POWERS;
      }

      rem = mmmMulmod(rem, powers[power], thePrime, _q);
#ifdef CHECK_CULLEN
      cullen = mmmSub(cullen, resC[power], thePrime);
#endif

#ifdef CHECK_WOODALL
      woodall = mmmSub(woodall, resW[power], thePrime);
#endif

      // At this point we have computed (1/b)^n (mod p).
      // If (1/b)^n (mod p) == n then we have a Woodall factor.
      // If (1/b)^n (mod p) == thePrime - n then we have a Cullen factor.

#ifdef CHECK_CULLEN
      if (rem == cullen)
         collect_factor(theN, +1, thePrime, factorCount, factors);
#endif

#ifdef CHECK_WOODALL
      if (rem == woodall)
         collect_factor(theN, -1, thePrime, factorCount, factors);
#endif

      prevN = theN;
      termIndex++;
   };
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

void collect_factor(int  n,
                    int  c,
                    long p,
                    volatile __global uint *factorCount,
                    __global long4 *factors)
{
   int old = atomic_inc(factorCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= MAX_FACTORS)
      return;

   factors[old].x = n;
   factors[old].y = c;
   factors[old].z = p;
}

long compute_inverse(long a, long p)
{
   long     m = p;
   long     y = 0, x = 1;
   long     q, t;

   while (a > 1)
   {
      // q is quotient
      q = a / m;
      t = m;

      // m is remainder now, process same as
      // Euclid's algo
      m = a % m;
      a = t;
      t = y;

      // Update y and x
      y = x - q * y;
      x = t;
   }

   // Make x positive
   if (x < 0)
      x += p;

   return x;
}
