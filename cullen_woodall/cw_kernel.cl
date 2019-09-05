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
                    volatile __global int *factorsCount,
                     __global long4 *factors);

long expmod(long base,
            long exp,
            ulong thePrime,
            ulong magicNumber,
            ulong magicShift);

long mulmod(long a,
            long b,
            ulong thePrime,
            ulong magicNumber,
            ulong magicShift);

long compute_inverse(long a, long modulus);

#define MAX_POWERS 50

__kernel void cw_kernel(__global const  long  *primes,
                        __global const ulong  *magicNumbers,
                        __global const ulong  *magicShifts,
                        __global const  int   *terms,
               volatile __global        int   *factorsCount,
                        __global        long4 *factors)
{
   int gid = get_global_id(0);

   int   termIndex;
   int   theN, prevN;
   int   power, idx;

   long  thePrime = primes[gid];
   ulong magicNumber = magicNumbers[gid];
   ulong magicShift = magicShifts[gid];
   long  rem = 0;
   ulong inverse;

   ulong powers[MAX_POWERS+1];
   
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
      
   theN = terms[0];
   
   rem = expmod(inverse, theN, thePrime, magicNumber, magicShift);
      
#ifdef CHECK_WOODALL
   if (rem == theN)
      collect_factor(theN, -1, thePrime, factorsCount, factors);
#endif

#ifdef CHECK_CULLEN
   if (rem == thePrime - theN)
      collect_factor(theN, +1, thePrime, factorsCount, factors);
#endif

   powers[1] = BASE;
      
   // Multiply successive terms by a (mod p)
   for (idx=2; idx<MAX_POWERS+1; idx++)
   {
      // If the bsae is odd, then all n must be even and thus the difference
      // between any two remaining n for the base must also be even.
      if ((idx > 2) && (BASE & 1))
      {
         if (idx & 1)
            continue;

         powers[idx] = mulmod(powers[idx-2], powers[2], thePrime, magicNumber, magicShift);
      }
      else 
      {
         powers[idx] = mulmod(powers[idx-1], powers[1], thePrime, magicNumber, magicShift);
      }
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
         rem = mulmod(rem, powers[MAX_POWERS], thePrime, magicNumber, magicShift);
         power -= MAX_POWERS;
      }
      
      rem = mulmod(rem, powers[power], thePrime, magicNumber, magicShift);
      
      // At this point we have computed (1/b)^n (mod p).
      // If (1/b)^n (mod p) == n then we have a Woodall factor.
      // If (1/b)^n (mod p) == thePrime - n then we have a Cullen factor.

#ifdef CHECK_WOODALL
      if (rem == theN)
         collect_factor(theN, -1, thePrime, factorsCount, factors);
#endif

#ifdef CHECK_CULLEN
      if (rem == thePrime - theN)
         collect_factor(theN, +1, thePrime, factorsCount, factors);
#endif

      prevN = theN;
      termIndex++;
   };
}

void collect_factor(int  n,
                    int  c,
                    long p,
                    volatile __global int *factorsCount,
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

long expmod(long base,
            long exp,
            ulong thePrime,
            ulong magicNumber,
            ulong magicShift)
{
   long x = base, y = 1;

   while (true)
   {
      if (exp & 1)
         y = mulmod(x, y, thePrime, magicNumber, magicShift);

      exp >>= 1;

      if (!exp)
         return y;

      x = mulmod(x, x, thePrime, magicNumber, magicShift);
   }
}

long mulmod(long a,
            long b,
            ulong thePrime,
            ulong magicNumber,
            ulong magicShift)
{
   ulong xa_low, xa_high;
   long xa_rem, xa_quot;
   ulong x196m, x196m1, x196m2, x196h;

   xa_low = a * b;
   xa_high = mul_hi(a, b);

   // xa_high | xa_low contains a 128-bit product of a*b

   x196m1 = mul_hi(xa_low, magicNumber);
   x196m2 = xa_high * magicNumber;
   x196h = mul_hi(xa_high, magicNumber);

   x196m = x196m1 + x196m2;
   if (x196m < x196m1) x196h++;

   xa_quot  = (x196m >> magicShift);
   xa_quot |= (x196h << (64 - magicShift));

   xa_rem = xa_low - (xa_quot * thePrime);
   if (xa_rem < 0) { xa_rem += thePrime; xa_quot -= 1; }

   return xa_rem;
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
