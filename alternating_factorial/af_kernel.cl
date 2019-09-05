/* afsieve.cl -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable
 
__kernel void af_kernel(__global const  long2 *primes,
                        __global const ulong2 *magicNumbers,
                        __global const ulong2 *magicShifts,
                        __global        long2 *rems,
                        __global        long2 *nm1rems,
                        __global        int   *n_range,
               volatile __global        long  *factors)
{
   int gid = get_global_id(0);

   long n_min = n_range[0];
   long n_start = n_min;
   long n_max = n_range[1];
   long maxMul;
   ulong2 n_mul;
   ulong2 quotient, xa_temp, xa_low, xa_high;
   ulong2 x196m1, x196m2, x196h, x196m;

   long2  pl = primes[gid];
   long2  xa_rem = rems[gid];
   long2  nm1_term = nm1rems[gid];
   ulong2 magicNumber = magicNumbers[gid];
   ulong2 magicShift = magicShifts[gid];

   if (pl.x == 0) return;
   if (pl.y == 0) pl.y = pl.x;

   maxMul = pl.y / n_max;

   n_mul.x = n_min;
   n_mul.y = n_min;

   do
   {
      if (xa_rem.x < maxMul && xa_rem.y < maxMul)
      {
         xa_rem.x = xa_rem.x * n_mul.x;
         xa_rem.y = xa_rem.y * n_mul.y;
      }
      else
      {
         xa_temp.x = xa_rem.x;
         xa_temp.y = xa_rem.y;
         xa_low = xa_temp * n_mul;
         xa_high = mul_hi(xa_temp, n_mul);

         // xa_high | xa_low contains a 128-bit product of a*b

         x196m1 = mul_hi(xa_low, magicNumber);
         x196m2 = xa_high * magicNumber;
         x196h = mul_hi(xa_high, magicNumber);

         x196m = x196m1 + x196m2;
         if (x196m.x < x196m1.x) x196h.x++;
         if (x196m.y < x196m1.y) x196h.y++;

         quotient = (x196m >> magicShift);
         quotient |= (x196h << (64 - magicShift));

         xa_rem.x = (long) (xa_low.x) - ((long) quotient.x * pl.x);
         xa_rem.y = (long) (xa_low.y) - ((long) quotient.y * pl.y);

         if (xa_rem.x < 0) xa_rem.x += pl.x;
         if (xa_rem.y < 0) xa_rem.y += pl.y;
      }

      // xa_rem = n! % p
      // nm1_term = af(n-1) % p      
      if (xa_rem.x == nm1_term.x) atom_min(&factors[n_min-n_start], pl.x);
      if (xa_rem.y == nm1_term.y) atom_min(&factors[n_min-n_start], pl.y);
         
      // Set nm1_term to af(n) % p
      if (xa_rem.x > nm1_term.x)
         nm1_term.x = xa_rem.x - nm1_term.x;
      else
         nm1_term.x = (xa_rem.x - nm1_term.x + pl.x);
      
      if (xa_rem.y > nm1_term.y)
         nm1_term.y = xa_rem.y - nm1_term.y;
      else
         nm1_term.y = (xa_rem.y - nm1_term.y + pl.y);
         
      n_mul.x += 1;
      n_mul.y += 1;

      n_min += 1;
   } while (n_min < n_max);
      
   rems[gid] = xa_rem;
   nm1rems[gid] = nm1_term;
}
