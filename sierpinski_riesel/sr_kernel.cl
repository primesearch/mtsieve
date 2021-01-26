/* cw_kernel.cl -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable

#define HASH_NOT_FOUND     UINT_MAX
#define HASH_MASK1         (1<<15)
#define HASH_MASK2         (HASH_MASK1-1)

ulong  invmod(ulong a, ulong p);

// return a (mod p), assuming -p < a < p
inline ulong  smod(long a, ulong p)
{
   ulong ua;
   if (a >= 0)
   {
      ua = (ulong) a;
      
      return (ua < p) ? ua : ua%p;
   }
   else
   {
      ua = (ulong) -a;
      
      return (ua < p) ? p-ua : p-(ua%p);
   }
}

inline ulong  lmod(long a, ulong p)
{
   if (a >= 0)
     return (a < p) ? a : a%p;
   else
     return (-a < p) ? p-(-a) : p-(-a%p);
}

// return a (mod p), assuming 0 <= a
inline ulong  umod(ulong a, ulong p)
{
   return (a < p) ? a : a%p;
}

ulong buildTables(ulong thePrime, ulong _q, ulong _r2, ulong _one, ulong *resBDCK64,
                        __global const ulong  *SEQ_K,
                        __global const  long  *SEQ_C,
                        __global const uint   *SUBSEQ_SEQ,
                        __global const uint   *SUBSEQ_Q);
ulong babySteps(ulong thePrime, ulong _q, ulong _r2, ulong _one, ushort *h_table, ushort *h_olist, ulong *h_BJ64);
void  hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64);
uint  hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64);

void collectFactor(ulong   p,
                   uint    ssIdx,
                   uint    n,
 volatile __global uint   *factorCount,
          __global ulong4 *factors);
                     
ulong mmmInvert(ulong p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);
   
__kernel void sr_kernel(__global const ulong  *primes,
                        __global const ulong  *SEQ_K,
                        __global const  long  *SEQ_C,
                        __global const uint   *SUBSEQ_SEQ,
                        __global const uint   *SUBSEQ_Q,
               volatile __global       uint   *factorCount,
                        __global       ulong4 *factors)
{
   int gid = get_global_id(0);
   
   // For our HashTable
   ushort  h_table[HASH_SIZE];
   ushort  h_olist[HASH_ELEMENTS];
   ulong   h_BJ64[HASH_ELEMENTS+1];
   ulong   resBDCK64[SUBSEQUENCES];

   ulong   thePrime = primes[gid];
   uint    idx, j, ssIdx;

   for (idx=0; idx<HASH_SIZE; idx++)
      h_table[idx] = 0;  
      
   for (idx=0; idx<HASH_ELEMENTS; idx++)
      h_BJ64[idx] = 0;  
   
   h_BJ64[HASH_ELEMENTS] = ULONG_MAX;   
   
   ulong _q = mmmInvert(thePrime);
   ulong _one = mmmOne(thePrime);
   ulong _r2 = mmmR2(thePrime, _q, _one);
   
   ulong resBM64 = buildTables(thePrime, _q, _r2, _one, resBDCK64, SEQ_K, SEQ_C, SUBSEQ_SEQ, SUBSEQ_Q);
   
   uint orderOfB = babySteps(thePrime, _q, _r2, _one, h_table, h_olist, h_BJ64);

   if (orderOfB > 0)
   {
      for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)
      {
         j = hashLookup(resBDCK64[ssIdx], h_table, h_olist, h_BJ64);
         
         while (j < SIEVE_RANGE)
         {
            collectFactor(thePrime, ssIdx, N_TERM(ssIdx, 0, j), factorCount, factors);
            
            j += orderOfB;
         }
      }
   }
   else
   {
      for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)
      {
         j = hashLookup(resBDCK64[ssIdx], h_table, h_olist, h_BJ64);
         
         if (j != HASH_NOT_FOUND)
            collectFactor(thePrime, ssIdx, N_TERM(ssIdx, 0, j), factorCount, factors);
      }
   }
   
   if (GIANT_STEPS < 2)
      return;
   
   resBM64 = mmmPowmod(resBM64, BABY_STEPS, thePrime, _q, _one);
   
   for (uint i=1; i<GIANT_STEPS; i++)
   {
      for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)
      {
         resBDCK64[ssIdx] = mmmMulmod(resBDCK64[ssIdx], resBM64, thePrime, _q);
         
         uint j = hashLookup(resBDCK64[ssIdx], h_table, h_olist, h_BJ64);
         
         if (j != HASH_NOT_FOUND)
            collectFactor(thePrime, ssIdx, N_TERM(ssIdx, i, j), factorCount, factors);
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
inline ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2)
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

ulong babySteps(ulong thePrime, ulong _q, ulong _r2, ulong _one, ushort *h_table, ushort *h_olist, ulong *h_BJ64)
{
   uint    j;
   ulong   resBase = mmmN(BASE, thePrime, _q, _r2);
   ulong   resBexpQ = mmmPowmod(resBase, BESTQ, thePrime, _q, _one);
   ulong   firstResBJ = mmmPowmod(resBexpQ, SIEVE_LOW, thePrime, _q, _one);
   ulong   resBJ = firstResBJ;

   for (j=0; j<BABY_STEPS; j++)
   {
      hashInsert(resBJ, j, h_table, h_olist, h_BJ64);
      
      resBJ = mmmMulmod(resBJ, resBexpQ, thePrime, _q);
      
      if (resBJ == firstResBJ)
         return j + 1;
   }
   
   return 0;
}

ulong buildTables(ulong thePrime, ulong _q, ulong _r2, ulong _one, ulong *resBDCK64,
                        __global const ulong  *SEQ_K,
                        __global const  long  *SEQ_C,
                        __global const uint   *SUBSEQ_SEQ,
                        __global const uint   *SUBSEQ_Q)
{
   uint  qIdx, seqIdx, ssIdx;
   ulong resBD64[BESTQ];
   ulong resCK64[SEQUENCES];
   
   ulong firstResBM64 = mmmN(invmod(BASE, thePrime), thePrime, _q, _r2);
   ulong resBM64 = firstResBM64;
   
   resBD64[0] = _one;
   
   for (qIdx=1; qIdx<BESTQ; qIdx++)
   {
      resBD64[qIdx] = resBM64;
      
      resBM64 = mmmMulmod(resBM64, firstResBM64, thePrime, _q);
   }
   
   for (seqIdx=0; seqIdx<SEQUENCES; seqIdx++)
   {
      ulong v1 = umod(SEQ_K[seqIdx], thePrime);
      
      ulong v2 = lmod(-SEQ_C[seqIdx], thePrime);
           
      ulong v3 = invmod(v1, thePrime);
      
      v2 = mmmN(v2, thePrime, _q, _r2);
      v3 = mmmN(v3, thePrime, _q, _r2);
      
      resCK64[seqIdx] = mmmMulmod(v2, v3, thePrime, _q);
   }

   // Compute -c/(k*b^d) (mod p) for each subsequence.
   for (ssIdx=0; ssIdx<SUBSEQUENCES; ssIdx++)
   {
      ulong ck = resCK64[SUBSEQ_SEQ[ssIdx]];
      ulong bd = resBD64[SUBSEQ_Q[ssIdx]];
      
      resBDCK64[ssIdx] = mmmMulmod(bd, ck, thePrime, _q);
   }
      
   return resBM64;
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
                   uint    ssIdx,
                   uint    n,
 volatile __global uint   *factorCount,
          __global ulong4 *factors)
{
   int old = atomic_inc(factorCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= MAX_FACTORS)
      return;
   
   factors[old].x = ssIdx;
   factors[old].y = n;
   factors[old].z = p;
}
