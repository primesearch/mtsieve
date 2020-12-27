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
                    volatile __global uint *factorsCount,
                     __global long4 *factors);

long compute_inverse(long a, long modulus);

ulong mmmInvert(ulong p);
ulong mmmOne(ulong _p);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p);
ulong mmmPowmod(ulong base, ulong exp, ulong _p, ulong _q);

#define MAX_POWERS 50

__kernel void cw_kernel(__global const  long  *primes,
                        __global const  int   *terms,
               volatile __global        uint  *factorsCount,
                        __global        long4 *factors)
{
   int gid = get_global_id(0);

   int   termIndex;
   int   theN, prevN;
   int   power, idx;

   long  thePrime = primes[gid];
   long  rem = 0;
   
#ifdef CHECK_CULLEN
   ulong cullen;
#endif

#ifdef CHECK_WOODALL
   ulong woodall;
#endif
   
   ulong inverse;

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

   // compute the inverse of b (mod p)
   inverse = compute_inverse(BASE, thePrime);
      
   ulong _q = mmmInvert(thePrime);
   
   theN = terms[0];
      
   rem = mmmPowmod(inverse, theN, thePrime, _q);

#ifdef CHECK_CULLEN
   cullen  = mmmN(thePrime - theN, thePrime);
   
   if (rem == cullen)
      collect_factor(theN, +1, thePrime, factorsCount, factors);
#endif

#ifdef CHECK_WOODALL
   woodall = mmmN(theN, thePrime);
   
   if (rem == woodall)
      collect_factor(theN, -1, thePrime, factorsCount, factors);
#endif

   powers[1] = mmmN(BASE, thePrime);
   powers[2] = mmmMulmod(powers[1], powers[1], thePrime, _q);

   resC[1] = mmmN(thePrime - 1, thePrime);
   resW[1] = mmmN(1, thePrime);
   
   resC[2] = mmmAdd(resC[1], resC[1], thePrime);
   resW[2] = mmmAdd(resW[1], resW[1], thePrime);
   
   if (BASE & 1)
   {
      // If the base is odd, then all n must be even and thus the difference
      // between any two remaining n for the base must also be even.
      for (idx=4; idx<MAX_POWERS+1; idx+=2)
      {
         powers[idx] = mmmMulmod(powers[idx-2], powers[2], thePrime, _q);
         resC[idx] = mmmAdd(resC[idx-2], resC[2], thePrime);
         resW[idx] = mmmAdd(resW[idx-2], resW[2], thePrime);
      }
   }
   else
   {
      powers[idx] = mmmMulmod(powers[idx-1], powers[1], thePrime, _q);
      resC[idx] = mmmAdd(resC[idx-1], resC[1], thePrime);
      resW[idx] = mmmAdd(resW[idx-1], resW[1], thePrime);
   }
      
   prevN = terms[0];
   
   // Note that the terms start at max N and decrease
   termIndex = 1;
   while (terms[termIndex] > 0)
   {
      theN = terms[termIndex];

      power = prevN - theN;

      // Our table isn't infinite in size, so we'll break down the power into smaller
      // chunks.  At worst this might cause an extra mulmod or two every once in a while.
      while (power > MAX_POWERS)
      {
         rem = mmmMulmod(rem, powers[MAX_POWERS], thePrime, _q);
         
         cullen = mmmSub(cullen, resC[MAX_POWERS], thePrime);
         woodall = mmmSub(woodall, resW[MAX_POWERS], thePrime);
         
         power -= MAX_POWERS;
      }
      
      rem = mmmMulmod(rem, powers[power], thePrime, _q);
      cullen = mmmSub(cullen, resC[power], thePrime);
      woodall = mmmSub(woodall, resW[power], thePrime);
      
      // At this point we have computed (1/b)^n (mod p).
      // If (1/b)^n (mod p) == n then we have a Woodall factor.
      // If (1/b)^n (mod p) == thePrime - n then we have a Cullen factor.

#ifdef CHECK_CULLEN
      if (rem == cullen)
         collect_factor(theN, +1, thePrime, factorsCount, factors);
#endif

#ifdef CHECK_WOODALL
      if (rem == woodall)
         collect_factor(theN, -1, thePrime, factorsCount, factors);
#endif

      prevN = theN;
      termIndex++;
   };
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

// Compute the residual of b ^ n (mod p)
ulong   mmmPowmod(ulong base, ulong exp, ulong _p, ulong _q)
{
   ulong x = mmmN(base, _p);
   ulong y = mmmN(1, _p);
   
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
                    volatile __global uint *factorsCount,
                     __global long4 *factors)
{
   int old = atomic_inc(factorsCount);

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
