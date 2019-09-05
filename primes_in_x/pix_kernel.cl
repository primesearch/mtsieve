/* pixsieve.cl -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

__kernel void pix_kernel(__global const  long2 *primes,
                         __global const ulong2 *magicNumbers,
                         __global const ulong2 *magicShifts,
                         __global        long2 *rems,
                         __global       uint   *digits,
                volatile __global        long  *factors)
{
   int gid = get_global_id(0);

   long maxMul;
   ulong2 multiplier;
   ulong2 quotient, xa_temp, xa_low, xa_high;
   ulong2 x196m1, x196m2, x196h, x196m;
   int  index = 0;

   long2  pl = primes[gid];
   long2  xa_rem = rems[gid];
   ulong2 magicNumber = magicNumbers[gid];
   ulong2 magicShift = magicShifts[gid];

   if (pl.x == 0) return;
   if (pl.y == 0) pl.y = pl.x;

   multiplier.x = digits[1];
   multiplier.y = digits[1];
   
   maxMul = pl.y / multiplier.y;

   index = 2;

   while (digits[index] < multiplier.x)
   {
      if (xa_rem.x < maxMul && xa_rem.y < maxMul)
      {
         // We don't need to do a mod in this case
         xa_rem.x = xa_rem.x * multiplier.x;
         xa_rem.y = xa_rem.y * multiplier.y;
      }
      else
      {
         xa_temp.x = xa_rem.x;
         xa_temp.y = xa_rem.y;
         xa_low = xa_temp * multiplier;
         xa_high = mul_hi(xa_temp, multiplier);

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

      xa_rem.x += digits[index];
      xa_rem.y += digits[index];
      
      if (xa_rem.x >= pl.x) xa_rem.x -= pl.x;
      if (xa_rem.y >= pl.y) xa_rem.y -= pl.y;
      
      if (xa_rem.x == 0) atom_min(&factors[index], pl.x);
      if (xa_rem.y == 0) atom_min(&factors[index], pl.y);

      index++;
   }
      
   rems[gid] = xa_rem;
}
