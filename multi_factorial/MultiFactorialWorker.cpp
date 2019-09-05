/* MultiFactorialWorker.cpp -- (C) Mark Rodenkirch, October 2012

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "MultiFactorialWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"

extern "C" int mfsieve(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);
extern "C" int multifactorial(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);

MultiFactorialWorker::MultiFactorialWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   uint32_t  maxN;
   
   ip_MultiFactorialApp = (MultiFactorialApp *) theApp;

   ib_Initialized = true;
   
   if (CpuSupportsAvx())
   {
      maxN = ip_MultiFactorialApp->GetMaxN();
      
      id_Terms = (double *) xmalloc((maxN + 1) * sizeof(double));
      
      for (uint32_t i=0; i<=maxN; i++)
         id_Terms[i] = (double) i;

      SetMiniChunkRange(maxN + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }   
}

void  MultiFactorialWorker::CleanUp(void)
{
   if (CpuSupportsAvx())
      xfree(id_Terms);
}

void  MultiFactorialWorker::TestMegaPrimeChunk(void)
{
   uint64_t  p[4], maxPrime = ip_App->GetMaxPrime();
   uint32_t  gotFactor, maxN = ip_MultiFactorialApp->GetMaxN();
   uint32_t  mf = ip_MultiFactorialApp->GetMultiFactorial();
   vector<uint64_t>::iterator it = iv_Primes.begin();

   while (it != iv_Primes.end())
   {
      p[0] = *it;
      it++;
      
      p[1] = *it;
      it++;
      
      p[2] = *it;
      it++;
      
      p[3] = *it;
      it++;

      for (uint32_t n=1; n<=mf; n++)
      {
         // If i is odd and mf is even, then i!mf is always odd, thus
         // n!mf+1 and n!mf-1 are always even.
         if (!(mf & 1) && (n & 1))
            continue;
         
         // The sieving routine requires that p > mf so that mf % p = mf.
         // As mf is small compared to p and as most small p will divide
         // at least one multi-factorial, this will have a negligable
         // affect on efficiency.
         if (p[0] > mf)
            gotFactor = mfsieve(n, mf, maxN, p);
         else 
            gotFactor = 1;
      
         if (gotFactor)
         {
            ExtractFactors(n, p[0]);
            ExtractFactors(n, p[1]);
            ExtractFactors(n, p[2]);
            ExtractFactors(n, p[3]);
         }
      }
      
      SetLargestPrimeTested(p[3], 4);
      
      if (p[3] >= maxPrime)
         break;
   }
}

void  MultiFactorialWorker::ExtractFactors(uint32_t start, uint64_t p)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.
   
   // We know that a factor was found, but the call to MultiFactorialWorker
   // doesn't tell us the term or the p for it, so we will use
   // the long way to find it.
   double   inverse, qd;
   uint32_t mf = ip_MultiFactorialApp->GetMultiFactorial();
   uint32_t maxN = ip_MultiFactorialApp->GetMaxN();
   uint32_t term, value;
   int64_t  q, r;
   long double pmax, vd;

   inverse = 1.0/p;
   r = 1;

   term = start;
   value = 1;
   
   pmax = (long double) (p + 1);
   vd = (long double) value;   
   
   while (term <= maxN)
   {
      vd *= (long double) term;
      value *= term;

      // If the current value is greater than the
      // the max, then the value cannot be prime.      
      if (vd > pmax)
         break;
      
      if (value == p+1) 
         ip_MultiFactorialApp->ReportPrime(p, term, -1);
         
      if (value == p-1) 
         ip_MultiFactorialApp->ReportPrime(p, term, +1);
            
      term += mf;
   }

   term = start;
   
   while (term <= maxN)
   {
      qd = ((double) r * (double) term);

      q = (int64_t) (qd * inverse);

      r = (r*term) - p*q;

      if (r < 0)
         r += p;
      else if (r >= (int64_t) p)
         r -= p;
      
      if (r == +1) 
         ip_MultiFactorialApp->ReportFactor(p, term, -1);
      
      if (r == (int64_t) p-1) 
         ip_MultiFactorialApp->ReportFactor(p, term, +1);
      
      term += mf;
   }
}

void  MultiFactorialWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   double __attribute__((aligned(32))) dps[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) reciprocals[AVX_ARRAY_SIZE];
   double __attribute__((aligned(32))) nextTerm[1];
   uint32_t term;
   uint32_t maxN = ip_MultiFactorialApp->GetMaxN();
   uint32_t mf = ip_MultiFactorialApp->GetMultiFactorial();
   
   // compute the inverse of b (mod p)
   for (int i=0; i<AVX_ARRAY_SIZE; i++)
      dps[i] = (double) miniPrimeChunk[i];
   
   avx_compute_reciprocal(dps, reciprocals);

   for (uint32_t n=1; n<=mf; n++)
   {
      nextTerm[0] = (double) n;
      
      avx_set_1a(nextTerm);
      
      for (term=n+mf; term<=maxN; term+=mf)
      {
         nextTerm[0] = id_Terms[term];
         avx_set_1b(nextTerm);
         avx_mulmod(dps, reciprocals);
         
         CheckAVXResult(miniPrimeChunk, dps, term);
      }
   }
}

void  MultiFactorialWorker::CheckAVXResult(uint64_t *ps, double *dps, uint32_t theN)
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
            VerifyAVXFactor(ps[idx], theN, -1);
   }
      
   // Only go further if one or more of the 16 primes yielded a factor for this n
   if (avx_neg_compare_1v(comparator, dps))
   {
      avx_get_16a(rems);
      
      for (idx=0; idx<AVX_ARRAY_SIZE; idx++)
         if (rems[idx] == dps[idx] - comparator[0])
            VerifyAVXFactor(ps[idx], theN, +1);
   }
}

void  MultiFactorialWorker::VerifyAVXFactor(uint64_t p, uint32_t theN, int32_t theC)
{
   // Yes, I know that this doesn't check for primes like ExtractFactors, but we know
   // all primorial primes < 2^128, so it isn't important that this doesn't catch them.
   uint64_t rem = 1;
   int32_t  n = (int32_t) theN;
   uint32_t mf = ip_MultiFactorialApp->GetMultiFactorial();
   bool     found = false;
   bool     modded = false;
   
   fpu_push_1divp(p);
   
   while (n > 1)
   {      
      if (rem * n > p + 1)
         modded = true;
      
      rem = fpu_mulmod(rem, n, p);
      
      n -= mf;
   }
   
   fpu_pop();
      
   if (rem == +1)
   {
      if (modded)
         ip_MultiFactorialApp->ReportFactor(p, theN, -1);
      else
         ip_MultiFactorialApp->ReportPrime(p, theN, -1);
   
      if (theC == -1)
         found = true;
   }
   
   if (p == rem + 1)
   {
      if (modded)
         ip_MultiFactorialApp->ReportFactor(p, theN, +1);
      else
         ip_MultiFactorialApp->ReportPrime(p, theN, +1);
   
      if (theC == +1)
         found = true;
   }   
 
   if (found)
      return;
   
   if (mf == 1)
      FatalError("The AVX routine found factor %" PRIu64" of %u!%+d, but it could not be validated", p, theN, theC);
   else
      FatalError("The AVX routine found factor %" PRIu64" of %u!%u%+d, but it could not be validated", p, theN, mf, theC);
}
