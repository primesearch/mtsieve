/* PrimorialWorker.cpp -- (C) Mark Rodenkirch, July 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "PrimorialWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"
extern "C" int primorial(const uint32_t *N, const uint64_t *P);

PrimorialWorker::PrimorialWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   uint32_t pIdx;
   
   ip_PrimorialApp = (PrimorialApp *) theApp;
   
   ip_PrimorialPrimes = ip_PrimorialApp->GetPrimorialPrimes(ii_NumberOfPrimorialPrimes);

   if (CpuSupportsAvx())
   {
      id_PrimorialPrimes = (double *) xmalloc((ii_NumberOfPrimorialPrimes + 1) * sizeof(double));
      
      for (pIdx=0; pIdx<ii_NumberOfPrimorialPrimes; pIdx++)
         id_PrimorialPrimes[pIdx] = (double) ip_PrimorialPrimes[pIdx];
      
      id_PrimorialPrimes[pIdx] = 0.0;
      SetMiniChunkRange(ip_PrimorialApp->GetMaxPrimorial() + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }
   
   ib_Initialized = true;
}

void  PrimorialWorker::CleanUp(void)
{
   if (id_PrimorialPrimes != NULL)
      xfree(id_PrimorialPrimes);
}

void  PrimorialWorker::TestMegaPrimeChunk(void)
{
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   vector<uint64_t>::iterator it = iv_Primes.begin();
      
   while (it != iv_Primes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;
   
      if (primorial(ip_PrimorialPrimes, ps))
      {
         ExtractFactors(ps[0]);
         ExtractFactors(ps[1]);
         ExtractFactors(ps[2]);
         ExtractFactors(ps[3]);
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] >= maxPrime)
         break;
   }
}

void  PrimorialWorker::ExtractFactors(uint64_t p)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.  This is okay since sieving past 2^52 is extremely
   // unlikely considering how slow sieving is for maxPrimorial.
   
   // We know that a factor was found, but the call to the ASM routine doesn't tell us the
   // term or the p for it, so we will use the long way to find it.
   double      inverse, qd;
   uint32_t    term, idx;
   int64_t     q, r;

   inverse = 1.0/p;
   r = FIRST_PRIMORIAL;
      
   for (idx=0; idx<ii_NumberOfPrimorialPrimes; idx++)
   {
      term = ip_PrimorialPrimes[idx];

      qd = ((double) r * (double) term);

      q = (int64_t) (qd * inverse);

      r = (r*term) - p*q;

      if (r < 0)
         r += p;
      else if (r >= (int64_t) p)
         r -= p;
      
      if (r == +1) 
         ip_PrimorialApp->ReportFactor(p, term, -1, false);
      
      if (r == (int64_t) p-1) 
         ip_PrimorialApp->ReportFactor(p, term, +1, false);
   }
}

void  PrimorialWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   double __attribute__((aligned(32))) dps[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) reciprocals[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) nextPrime[1];
   uint32_t idx;
   
   // compute the inverse of b (mod p)
   for (int i=0; i<AVX_ARRAY_SIZE; i++)
      dps[i] = (double) miniPrimeChunk[i];
   
   avx_compute_reciprocal(dps, reciprocals);
   
   nextPrime[0] = (double) FIRST_PRIMORIAL;
   
   avx_set_1a(nextPrime);

   for (idx=0; idx<ii_NumberOfPrimorialPrimes; idx++)
   {
      nextPrime[0] = id_PrimorialPrimes[idx];
      
      avx_set_1b(nextPrime);
      avx_mulmod(dps, reciprocals);
      
      CheckAVXResult(miniPrimeChunk, dps, ip_PrimorialPrimes[idx]);
   }
}

void  PrimorialWorker::CheckAVXResult(uint64_t *ps, double *dps, uint32_t primorial)
{
   uint32_t idx;
   double __attribute__((aligned(32))) comparator[1];
   double __attribute__((aligned(32))) rems[AVX_ARRAY_SIZE];
   
   comparator[0] = 1.0;
   
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (avx_pos_compare_1v(comparator) > 0)
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == comparator[0])
            ip_PrimorialApp->ReportFactor(ps[idx], primorial, -1, true);
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (avx_neg_compare_1v(comparator, dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - comparator[0])
            ip_PrimorialApp->ReportFactor(ps[idx], primorial, +1, true);
   }
}
