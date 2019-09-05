/* xyyxsieve.cl -- (C) Mark Rodenkirch, June 2014

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */
 
#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

#define MAX_POWERS 50

void collect_factor(int  x,
                    int  y,
                    int  c,
                    long p,
                    volatile __global int *counter,
                     __global long4 *factors);

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

__kernel void xyyx_kernel(__global const  long  *primes,
                          __global const ulong  *magicNumbers,
                          __global const ulong  *magicShifts,
                          __global const  int   *terms,
                volatile  __global        long  *minusFactors,
                volatile  __global        long  *plusFactors)
{
   int gid = get_global_id(0);

   int   termIndex;
   int   currentX;
   int   currentY, previousY;
   int   pIndex;

   long  thePrime = primes[gid];
   ulong magicNumber = magicNumbers[gid];
   ulong magicShift = magicShifts[gid];
   long  xPowY, xPowN;
   long  yPowX;

   ulong powers[MAX_POWERS];

   if (thePrime == 0)
      return;

   termIndex = 0;

   while (terms[termIndex] > 0)
   {
      currentX = terms[termIndex];
      termIndex++;
      
      previousY = 0;

      // If we have at least one y for this x
      if (terms[termIndex] > 0)
      {
         currentY = terms[termIndex];

         xPowY = mulmod(currentX, currentX, thePrime, magicNumber, magicShift);
                  
         powers[0] = 1;
         powers[1] = xPowY;
         
         // Multiply successive terms by x^2 (mod p)
         for (pIndex=2; pIndex<=MAX_POWERS; pIndex++)
             powers[pIndex] = mulmod(powers[pIndex-1], xPowY, thePrime, magicNumber, magicShift);
       
         xPowY = powmod(currentX, currentY, thePrime, magicNumber, magicShift);

         previousY = currentY;
  
         while (terms[termIndex] > 0)
         {
            currentY = terms[termIndex];
      
            // If x is even then y must be odd and if x is odd
            // then y must be even so y - prevY is always even
            pIndex = (currentY - previousY) >> 1;
            
            if (pIndex >= MAX_POWERS)
               xPowN = powmod(currentX, currentY - previousY, thePrime, magicNumber, magicShift);
            else
               xPowN = powers[pIndex];

            // We have x^prevY (mod p).
            // Now compute x^y (mod p) as (x^prevY * x^(y-prevY) (mod p)
            xPowY = mulmod(xPowY, xPowN, thePrime, magicNumber, magicShift);
            
            yPowX = powmod(currentY, currentX, thePrime, magicNumber, magicShift);

            // x^y - y^x = 0(mod p)
            if (xPowY == yPowX) atom_min(&minusFactors[termIndex], thePrime);

            // x^y + y^x == 0(mod p)
            if (xPowY + yPowX == thePrime) atom_min(&plusFactors[termIndex], thePrime);
               
            termIndex++;
            previousY = currentY;
         }
      }
      
      // Get the next x
      termIndex++;
   };
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
