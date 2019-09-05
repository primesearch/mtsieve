
/* magic.cl -- (C) Mark Rodenkirch, December 2013

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

void compute_magic(long theDivisor,
                   ulong *magicNumber,
                   ulong *magicShift);
                   
long powmod(long base,
            long exp,
            ulong thePrime,
            ulong magicNumber,
            ulong magicShift);
            
long mulmod(long a,
            long b,
            ulong thePrime,
            ulong magicNumber,
            ulong magicShift);

__kernel void gfn_kernel(__global const long *primes,
                         __global      ulong *ks)
{
   int    gid = get_global_id(0);

   long   p = primes[gid];
   long   k, rem;
   ulong  magicNumber, magicShift;
   ulong  mN, mS;
   int    n, index;

   compute_magic(p, &magicNumber, &magicShift);
   
   k = (1 + p) >> 1;

   rem = powmod(k, N_MIN, p, magicNumber, magicShift);
   
   index = gid * N_COUNT;
   
   k = p - rem;

   // Note that p > K_MAX is a requirement when using this code
   for (n=0; n<N_COUNT; n++)
   {    
      if (k & 1)
      {
         if (k <= K_MAX && k >= K_MIN)
            ks[index + n] = k;

         k += p;
      }
      else
         ks[index + n] = 0;
      
      // For the next n, divide k by 2
      // Note that k*2^n+1 = (k/2)*2^(n+1)+1

      k >>= 1;
   }
}

void compute_magic(long theDivisor,
                   ulong *magicNumber,
                   ulong *magicShift)
{
   ulong two63 = 0x8000000000000000;

   ulong d = theDivisor;
   ulong t = two63;
   ulong anc = t - 1 - t%d;    // Absolute value of nc.
   ulong p = 63;               // Init p.
   ulong q1 = two63/anc;       // Init q1 = 2**p/|nc|.
   ulong r1 = two63 - q1*anc;  // Init r1 = rem(2**p, |nc|).
   ulong q2 = two63/d;         // Init q2 = 2**p/|d|.
   ulong r2 = two63- q2*d;     // Init r2 = rem(2**p, |d|).
   ulong delta;

   do {
      p = p + 1;
      q1 = 2*q1;               // Update q1 = 2**p/|nc|.
      r1 = 2*r1;               // Update r1 = rem(2**p, |nc|.
      if (r1 >= anc) {         // Must be an unsigned comparison
         q1 = q1 + 1; 
         r1 = r1 - anc;
      }
      q2 = 2*q2;               // Update q2 = 2**p/|d|.
      r2 = 2*r2;               // Update r2 = rem(2**p, |d|.
      if (r2 >= d) {           // Must be an unsigned comparison
         q2 = q2 + 1; 
         r2 = r2 - d;
      }
      delta = d - r2;
  } while (q1 < delta || (q1 == delta && r1 == 0));

   *magicNumber = (q2 + 1);
   *magicShift = (p - 64);
}

long powmod(long base,
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
