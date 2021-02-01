/* cw_kernel.cl -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

#define SP_MIXED     0
#define SP_EVEN      1
#define SP_ODD       2
#define SP_NO_PARITY 999

#define N_TERM(q, i, j)       ((SIEVE_LOW + (j) + (i)*bSteps)*BESTQ + q)
#define CSS_INDEX(x, y, z)    (((((x) * PRL_COUNT) + (y)) * POWER_RESIDUE_LCM) + (z))

#define HASH_NOT_FOUND        UINT_MAX
#define HASH_MASK1            (1<<15)
#define HASH_MASK2            (HASH_MASK1-1)

ushort getParity(ulong thePrime);
ulong  invmod(ulong a, ulong p);
short  legendre(long a, ulong p);

uint  setupDiscreteLog(ulong thePrime, ulong _q, ulong _one,
                       ulong resBase, ulong resInvBase, ulong resNegCK, ushort parity,
                       __global const short  *divisorShifts,
                       __global const uint   *prlIndices);

ulong buildLookupsAndClimbLadder(ulong thePrime, ulong _q, ulong _one,
                                  ulong resBase, ulong resNegCK, ulong *resBD,
                                  uint qIndex, uint ladderIndex,
                                  __global const ushort *qs,
                                  __global const ushort *ladders);

ulong babySteps(ulong thePrime, ulong _q, ulong _one, ulong resInvBase, uint bSteps,
                ushort *h_table, ushort *h_olist, ulong *h_BJ64);
                
void  hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64);
uint  hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64);

void collectFactor(ulong   p,
                   uint    n,
 volatile __global uint   *factorCount,
          __global ulong2 *factors);

ulong mmmInvert(ulong p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);
   
__kernel void cisonesingle_kernel(__global const ulong  *primes,
                                  __global const uint   *babyStepsArray,
                                  __global const uint   *giantStepsArray,
                                  __global const short  *divisorShifts,
                                  __global const uint   *prlIndices,
                                  __global const uint   *qIndices,
                                  __global const ushort *qs,
                                  __global const uint   *ladderIndices,
                                  __global const ushort *ladders,
                         volatile __global       uint   *factorCount,
                                  __global       ulong2 *factors)
{
   int    gid = get_global_id(0);
   ulong  thePrime = primes[gid];
   ulong  negCK;
   ushort parity;
   
   parity = getParity(thePrime);
      
   if (parity == SP_NO_PARITY)
      return;
   
   // neg_ck <-- -k/c (mod p) == -ck (mod p)
   if (thePrime < SEQ_K)
      negCK = SEQ_K - (SEQ_K / thePrime);
   else 
      negCK = SEQ_K;
  
   if (SEQ_C > 0)
      negCK = thePrime - negCK;
   
   ulong   _q = mmmInvert(thePrime);
   ulong   _one = mmmOne(thePrime);
   ulong   _r2 = mmmR2(thePrime, _q, _one);
   
   // compute 1/base (mod p)
   ulong   invBase = invmod(BASE, thePrime);

   ulong   resBase = mmmNToRes(BASE, thePrime, _q, _r2);
   ulong   resInvBase = mmmNToRes(invBase, thePrime, _q, _r2);
   ulong   resNegCK = mmmNToRes(negCK, thePrime, _q, _r2);
   ulong   resBexpQ;
   ushort  h_table[HASH_SIZE];
   ushort  h_olist[HASH_ELEMENTS];
   ulong   h_BJ64[HASH_ELEMENTS+1];
   ulong   resBD[SUBSEQUENCE_COUNT];
   uint    i, j, k, ssCount;
   uint    idx, qIdx, ladderIdx;

   ushort cssIndex = setupDiscreteLog(thePrime, _q, _one, resBase, resInvBase, resNegCK, parity, divisorShifts, prlIndices);
   
   qIdx = qIndices[cssIndex];

   // If no qs for this p, then no factors, so return
   if (qIdx == 0)
      return;
   
   ladderIdx = ladderIndices[cssIndex];
      
   // -ckb^d is an r-th power residue for at least one term (k*b^d)*(b^Q)^(n/Q)+c of this subsequence
   resBexpQ = buildLookupsAndClimbLadder(thePrime, _q, _one, resBase, resNegCK, resBD, qIdx, ladderIdx, qs, ladders);
   
   ssCount = qs[qIdx];

   for (idx=0; idx<HASH_SIZE; idx++)
      h_table[idx] = 0;  
      
   for (idx=0; idx<HASH_ELEMENTS; idx++)
      h_BJ64[idx] = 0;  
   
   h_BJ64[HASH_ELEMENTS] = ULONG_MAX;   
   
   uint bSteps = babyStepsArray[ssCount - 1];
   uint gSteps = giantStepsArray[ssCount - 1];

   ushort orderOfB = babySteps(thePrime, _q, _one, resInvBase, bSteps, h_table, h_olist, h_BJ64);

   if (orderOfB > 0)
   {
      // If orderOfB > 0, then this is all the information we need to
      // determine every solution for this p, so no giant steps are neede
      for (k=0; k<ssCount; k++)
      {
         j = hashLookup(resBD[k], h_table, h_olist, h_BJ64);

         while (j < bSteps * gSteps)
         {
            collectFactor(thePrime, N_TERM(qs[qIdx+1+k], 0, j), factorCount, factors);
             
            j += orderOfB;
         }
      }
   }
   else
   {
      // First giant step
      for (k=0; k<ssCount; k++)
      {
         j = hashLookup(resBD[k], h_table, h_olist, h_BJ64);

         if (j != HASH_NOT_FOUND)
            collectFactor(thePrime, N_TERM(qs[qIdx+1+k], 0, j), factorCount, factors);
      }

      // Remaining giant steps
      if (gSteps > 1)
      {
         ulong resBQM = mmmPowmod(resBexpQ, bSteps, thePrime, _q, _one);
         
         for (i=1; i<gSteps; i++)
         {
            for (k=0; k<ssCount; k++)
            {
               resBD[k] = mmmMulmod(resBD[k], resBQM, thePrime, _q);
               
               j = hashLookup(resBD[k], h_table, h_olist, h_BJ64);

               if (j != HASH_NOT_FOUND)
                  collectFactor(thePrime, N_TERM(qs[qIdx+1+k], i, j), factorCount, factors);
            }
         }
      }
   }
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
inline ulong mmmOne(ulong _p)
{
   return ((-_p) % _p);
}

// Compute the residual of 2^64 (mod p)
inline ulong mmmR2(ulong _p, ulong _q, ulong _one)
{
	ulong t = mmmAdd(_one, _one, _p);
   
   t = mmmAdd(t, t, _p);   // 4
	for (size_t i=0; i<5; i++)
      t = mmmMulmod(t, t, _p, _q);   // 4^{2^5} = 2^64
      
	return t;
}

inline ulong mmmAdd(ulong a, ulong b, ulong _p)
{
   ulong c = (a >= _p - b) ? _p : 0;
   return a + b - c;
}

inline ulong mmmSub(ulong a, ulong b, ulong _p)
{
   ulong c = (a < b) ? _p : 0;
   return a - b + c;
}

// Compute the residual of n (mod p)
inline ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2)
{
   return mmmMulmod(n, _r2, _p, _q);
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

// Compute the residual of b ^ n (mod p)
ulong   mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one)
{
   ulong x = resbase;
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

ushort  getParity(ulong thePrime)
{
   // Mixed parity sequences
   if (SEQ_PARITY == SP_MIXED)
   {
      short qr_m1, qr_p1;
      
     //if (ib_UseLegendreTables)
     //{
     //   uint32_t qr_mod = (thePrime/2) % ip_Legendre->mod;
     //            
     //   qr_m1 = (ip_Legendre->dualParityMapM1[L_BYTE(qr_mod)] & L_BIT(qr_mod));
     //   qr_p1 = (ip_Legendre->dualParityMapP1[L_BYTE(qr_mod)] & L_BIT(qr_mod));
     //}
     //else
      {
         short sym = legendre(KC_CORE, thePrime);
         
         qr_m1 = (sym == 1);
         qr_p1 = (sym == legendre(BASE, thePrime));
      }
      
      if (qr_m1)
         return (qr_p1 ? SP_MIXED : SP_EVEN);

      return (qr_p1 ? SP_ODD : SP_NO_PARITY);
   }
   
   // Single parity sequences
   short qr;
   
   //if (ib_UseLegendreTables)
   //{
   //   uint32_t qr_mod = (p/2) % ip_Legendre->mod;
   //      
   //   qr = (ip_Legendre->oneParityMap[L_BYTE(qr_mod)] & L_BIT(qr_mod));
   //}
   //else
   {
      short sym = legendre(KC_CORE, thePrime);
      
      if (SEQ_PARITY == SP_EVEN)
         qr = (sym == 1);
      else
         qr = (sym == legendre(BASE, thePrime));
   }
   
   if (qr)
      return SEQ_PARITY;
      
   return SP_NO_PARITY;
}

uint setupDiscreteLog(ulong thePrime, ulong _q, ulong _one,
                      ulong resBase, ulong resInvBase, ulong resNegCK, ushort parity,
                      __global const short  *divisorShifts,
                      __global const uint   *prlIndices)
{
   ulong pShift;
   uint  idx;
   uint  h, r;
   short shift;
   ulong resX[POWER_RESIDUE_LCM + 1];
   
   idx = (thePrime/2) % (POWER_RESIDUE_LCM/2);
   shift = divisorShifts[idx];
      
   if (shift == 0)
   {
      r = prlIndices[1];
      h = 0;
   
      return CSS_INDEX(parity, r, h);
   }
   
   if (shift > 0)
   {
      // p = 1 (mod s), where s is not a power of 2. Check for r-th power
      // residues for each prime power divisor r of s.
      pShift = thePrime / shift;
   }
   else
   {
      // p = 1 (mod 2^s), where s > 1. Check for r-th power residues for each divisor
      // r of s. We handle this case seperately to avoid computing p/s using plain division.
      pShift = thePrime >> (-shift);
      shift = 1 << (-shift);
   }
     
   resX[0] = _one;
   resX[1] = mmmPowmod(resInvBase, pShift, thePrime, _q, _one);
  
   for (r=1; resX[r] != resX[0]; r++)
      resX[r+1] = mmmMulmod(resX[r], resX[1], thePrime, _q);

   if (shift % r != 0)
      return 0;
    
   // resX[r] <- (-ck)^((p-1)/shift)
   resX[r] = mmmPowmod(resNegCK, pShift, thePrime, _q, _one);

   // Find h such that resX[h] = resX[r], i.e. (-ckb^h)^((p-1)/r)=1 (mod p), or h=r if not found
   for (h=0; resX[r] != resX[h]; h++)
      ;

   // If no h was found, then there is nothing further to do.
   if (h == r)
      return 0;

   r = prlIndices[r];
   
   return CSS_INDEX(parity, r, h);
}

ulong buildLookupsAndClimbLadder(ulong thePrime, ulong _q,  ulong _one,
                                  ulong resBase, ulong resNegCK, ulong *resBD,
                                  uint qIndex, uint ladderIndex,
                                  __global const ushort *qs,
                                  __global const ushort *ladders)
{
   uint   i, j, idx, lLen, qLen;
   ulong  resX[POWER_RESIDUE_LCM + 1];

   lLen = ladders[ladderIndex];
   
   // Precompute b^d (mod p) for 0 <= d <= Q, as necessary
   resX[0] = _one;
   resX[1] = resBase;
   resX[2] = mmmMulmod(resX[1], resX[1], thePrime, _q);
   
   i = 2;
   for (j=0; j<lLen; j++)
   {
      idx = ladders[ladderIndex+j+1];
      
      resX[i+idx] = mmmMulmod(resX[i], resX[idx], thePrime, _q);
      
      i += idx;
   }

   qLen = qs[qIndex];
   
   for (j=0; j<qLen; j++)
   {
      idx = qs[qIndex+j+1];
      
      resBD[j] = mmmMulmod(resX[idx], resNegCK, thePrime, _q);
   }

   return resX[BESTQ];
}

ulong babySteps(ulong thePrime, ulong _q, ulong _one,
                ulong resInvBase, uint bSteps,
                ushort *h_table, ushort *h_olist, ulong *h_BJ64)
{
   ulong    resInvBaseExpQ;
   ulong    resBJ;
   ulong    firstResBJ;
   uint     j;
    
   // b <- inv_b^Q (mod p)
   resInvBaseExpQ = mmmPowmod(resInvBase, BESTQ, thePrime, _q, _one);
  
   firstResBJ = resBJ = mmmPowmod(resInvBaseExpQ, SIEVE_LOW, thePrime, _q, _one);

   for (j=0; j<bSteps; j++)
   { 
      hashInsert(resBJ, j, h_table, h_olist, h_BJ64);
      
      resBJ = mmmMulmod(resBJ, resInvBaseExpQ, thePrime, _q);
            
      if (resBJ == firstResBJ)
         return j + 1;
   }
   
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

short  legendre(long a, ulong p)
{
   ulong x, y, t;
   short sign;

   if (a < 0)
   {
      a = -a;
      sign = (p % 4 == 1) ? 1 : -1;
   }
   else
      sign = 1;

   for (y = a; y % 2 == 0; y /= 2)
      if (p % 8 == 3 || p % 8 == 5)
         sign = -sign;

   if (p % 4 == 3 && y % 4 == 3)
      sign = -sign;

   for (x = p % y; x > 0; x %= y)
   {
      for ( ; x % 2 == 0; x /= 2)
         if (y % 8 == 3 || y % 8 == 5)
            sign = -sign;

      t = x, x = y, y = t;

      if (x % 4 == 3 && y % 4 == 3)
         sign = -sign;
   }

   return sign;
}

void hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64)
{
   uint slot;

   h_BJ64[j] = bj;
   slot = bj & (HASH_SIZE - 1);
   if (h_table[slot] == HASH_ELEMENTS)
      h_table[slot] = j;
   else
   {
      h_olist[j] = (h_table[slot] ^ HASH_MASK1);
      h_table[slot] = (j | HASH_MASK1);
   }
}

uint hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64)
{
   uint   slot;
   ushort elt;

   slot = bj & (HASH_SIZE - 1);
   elt = h_table[slot];

   if (h_BJ64[elt & HASH_MASK2] == bj)
      return elt & HASH_MASK2;
   
   if ((elt & HASH_MASK1) == 0)
      return HASH_NOT_FOUND;

   elt &= HASH_MASK2;
   do
   {
      if (h_BJ64[h_olist[elt] & HASH_MASK2] == bj)
         return h_olist[elt] & HASH_MASK2;
      
      elt = h_olist[elt];
   } while ((elt & HASH_MASK1) == 0);

   return HASH_NOT_FOUND;
}

void collectFactor(ulong   p,
                   uint    n,
 volatile __global uint   *factorCount,
          __global ulong2 *factors)
{
   int old = atomic_inc(factorCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= MAX_FACTORS)
      return;
   
   factors[old].x = n;
   factors[old].y = p;
}
