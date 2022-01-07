/* mf_kernel.cl -- (C) Mark Rodenkirch, September 2016

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
 */

#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

ulong  invmod(ulong a, ulong p);

void collect_factor(uint    n,
                    ulong   p,
  volatile __global uint   *factorsCount,
           __global ulong2 *factors);

ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one);
ulong mmmResToN(ulong res, ulong _p, ulong _q);

__kernel void sm_kernel(__global const  ulong  *primes,
                        __global const  uint   *terms,
               volatile __global        uint   *factorsCount,
                        __global        ulong2 *factors)
{
   int    gid = get_global_id(0);
   
   ulong  n;
   ulong  thePrime = primes[gid];

   if (thePrime == 0) return;
   
   ulong _q = mmmInvert(thePrime);
   ulong _one = mmmOne(thePrime);
   ulong _r2 = mmmR2(thePrime, _q, _one);

   ulong invmod2 = invmod(12321, thePrime);
   ulong invmod3 = invmod(1234321, thePrime);
   ulong invmod4 = invmod(123454321, thePrime);
      
   ulong resTenE1 = mmmNtoRes(10, thePrime, _q, _r2);
      
   ulong tempMul1 = mmmNtoRes(15208068915062105958, thePrime, _q, _r2);
   ulong tempSub1 = mmmNtoRes(11211123422, thePrime, _q, _r2);
   
   ulong tempSub2 = mmmNtoRes(1109890222, thePrime, _q, _r2);
   ulong resInvmod2 = mmmNtoRes(invmod2, thePrime, _q, _r2);
   
   ulong tempMul3 = mmmNtoRes(123454321, thePrime, _q, _r2);
   ulong tempSub3 = mmmNtoRes(1110988902222, thePrime, _q, _r2);
   ulong resInvmod3 = mmmNtoRes(invmod3, thePrime, _q, _r2);
   
   ulong tempMul4 = mmmNtoRes(12345654321, thePrime, _q, _r2);
   ulong tempSub4 = mmmNtoRes(1111098889022222, thePrime, _q, _r2);
   ulong resInvmod4 = mmmNtoRes(invmod4, thePrime, _q, _r2);
   
   //     t = powmod(10,179,f);
   //     C = mulmod(t,15208068915062105958%f,f);
   //     submod(C,11211123422);
   
   ulong resTemp = mmmPowmod(resTenE1, 179, thePrime, _q, _one);
   ulong resC = mmmMulmod(resTemp, tempMul1, thePrime, _q);
   resC = mmmSub(resC, tempSub1, thePrime);

   //     t = powmod(10,2699,f);
   //     C = mulmod(C, t, f);
   //     submod(C,1109890222);
   //     divmod(C, 12321);

   resTemp = mmmPowmod(resTenE1, 2699, thePrime, _q, _one);
   resC = mmmMulmod(resC, resTemp, thePrime, _q);
   resC = mmmSub(resC, tempSub2, thePrime);
   resC = mmmMulmod(resC, resInvmod2, thePrime, _q);

   //     t = powmod(10,35999,f);
   //     C = mulmod(C, t, f);
   //     C = mulmod(C, 123454321, f);
   //     submod(C,1110988902222);
   //     divmod(C, 1234321);
   
   resTemp = mmmPowmod(resTenE1, 35999, thePrime, _q, _one);
   resC = mmmMulmod(resC, resTemp, thePrime, _q);
   resC = mmmMulmod(resC, tempMul3, thePrime, _q);
   resC = mmmSub(resC, tempSub3, thePrime);
   resC = mmmMulmod(resC, resInvmod3, thePrime, _q);
   resC = mmmMulmod(resC, tempMul4, thePrime, _q);
   
   //      C = mulmod(C,12345654321 % f, f);
   //      t = powmod(10,449999,f);
   //      C = mulmod(C, t, f);
   //      submod(C,1111098889022222);
   //      divmod(C, 123454321);
   
   resTemp = mmmPowmod(resTenE1, 449999, thePrime, _q, _one);
   resC = mmmMulmod(resC, resTemp, thePrime, _q);
   resC = mmmSub(resC, tempSub4, thePrime);
   resC = mmmMulmod(resC, resInvmod4, thePrime, _q);
   
   //      t = powmod(10,6*v[0]-599989,f);
   //      C = mulmod(C, t, f);
   //      M = (v[0]*999999+1000000)%f;
   //      r = 1000000000000000000%f;
   
   uint   exp = 6*terms[0] - 599989;
   ulong  m = ((ulong) terms[0] * 999999) + 1000000;
   
   resTemp = mmmPowmod(resTenE1, exp, thePrime, _q, _one);
   resC = mmmMulmod(resC, resTemp, thePrime, _q);
   
   ulong resM = mmmNtoRes(m, thePrime, _q, _r2);
   
   ulong res9s = mmmNtoRes(999999, thePrime, _q, _r2);
   ulong resR = mmmNtoRes(1000000000000000000, thePrime, _q, _r2);
   ulong resT[6];

   //      T[1] = mulmod(r, r, f);
   //      T[2] = mulmod(T[1], T[1], f);
   //      T[3] = mulmod(T[1], T[2], f);
   //      T[4] = mulmod(T[2], T[2], f);
   //      T[5] = mulmod(T[2], T[3], f);

   resT[1] = mmmMulmod(resR, resR, thePrime, _q);
   resT[2] = mmmMulmod(resT[1], resT[1], thePrime, _q);
   resT[3] = mmmMulmod(resT[1], resT[2], thePrime, _q);
   resT[4] = mmmMulmod(resT[2], resT[2], thePrime, _q);
   resT[5] = mmmMulmod(resT[2], resT[3], thePrime, _q);
   
   // This is only for the first term    
      if (resC == resM)
         collect_factor(terms[0], thePrime, factorsCount, factors);

   // This is for the remaining terms
   uint idx = 1;
   while (terms[idx] > 0)
   {         
      uint dn = terms[idx] - terms[idx-1];
      
      //    if(dn<=30)
      //      C = mulmod(C, T[dn/6], f);
      //    else {
      //      t = powmod(r, dn/3, f);
      //      C = mulmod(C, t, f);
      //    }
      
      if (dn <= 30)
         resC = mmmMulmod(resC, resT[dn/6], thePrime, _q);
      else
      {
         resTemp = mmmPowmod(resT[1], dn/6, thePrime, _q, _one);
         resC = mmmMulmod(resC, resTemp, thePrime, _q);
      }
      
      // M = (M + 999999*dn) % f;
      
      resTemp = mmmNtoRes(dn, thePrime, _q, _r2);
      resTemp = mmmMulmod(resTemp, res9s, thePrime, _q);
      resM = mmmAdd(resM, resTemp, thePrime);
      
      if (resC == resM)
         collect_factor(terms[idx], thePrime, factorsCount, factors);
         
      idx++;
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
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one)
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

ulong  invmod(ulong a, ulong p)
{
   ulong ps1, ps2, q, r, t, dividend, divisor;
   uint parity;

   if (a < 3)
      return (a < 2) ? a : (p+1)/2;

   q = p / a;
   r = p % a;

   dividend = a;
   divisor = r;
   ps1 = q;
   ps2 = 1;
   parity = 0;

   while (divisor > 1)
   {
      r = dividend - divisor;
      t = r - divisor;
      if (r >= divisor) {
         q += ps1; r = t; t -= divisor;
         if (r >= divisor) {
            q += ps1; r = t; t -= divisor;
            if (r >= divisor) {
               q += ps1; r = t; t -= divisor;
               if (r >= divisor) {
                  q += ps1; r = t; t -= divisor;
                  if (r >= divisor) {
                     q += ps1; r = t; t -= divisor;
                     if (r >= divisor) {
                        q += ps1; r = t; t -= divisor;
                        if (r >= divisor) {
                           q += ps1; r = t; t -= divisor;
                           if (r >= divisor) {
                              q += ps1; r = t;
                              if (r >= divisor) {
                                 q = dividend / divisor;
                                 r = dividend % divisor;
                                 q *= ps1;
                              } } } } } } } } }
      q += ps2;
      parity = ~parity;
      dividend = divisor;
      divisor = r;
      ps2 = ps1;
      ps1 = q;
   }

   return (parity) ? ps1 : p - ps1;
}

void collect_factor(uint    n,
                    ulong   p,
  volatile __global uint   *factorsCount,
           __global ulong2 *factors)
{

   uint old = atomic_inc(factorsCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= D_MAX_FACTORS)
      return;
   
   factors[old].x = n;
   factors[old].y = p;
}