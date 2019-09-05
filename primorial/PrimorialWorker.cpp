/* PrimorialWorker.cpp -- (C) Mark Rodenkirch, July 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "PrimorialWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/avx-asm-x86.h"
#include "../sieve/primesieve.hpp"

extern "C" int primorial(const uint32_t *N, const uint64_t *P);

PrimorialWorker::PrimorialWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   uint32_t primesToTest, idx;
   vector<uint64_t>  primes;
   vector<uint64_t>::iterator it;
   
   ip_PrimorialApp = (PrimorialApp *) theApp;
   
   primesToTest = primesieve::count_primes(0, ip_PrimorialApp->GetMaxN() + 1);
   
   ii_Primes = (uint32_t *) xmalloc((primesToTest + 1) * sizeof(uint32_t));
   
   if (CpuSupportsAvx())
      id_Primes = (double *) xmalloc((primesToTest + 1) * sizeof(double));

   primes.clear();
   
   // Generate primes for this worker
   primesieve::generate_n_primes(primesToTest, 0, &primes);
   
   it = primes.begin();
   idx = 0;
   
   while (it != primes.end())
   {
      ii_Primes[idx] = *it;
      
      if (CpuSupportsAvx())
         id_Primes[idx] = (double) *it;
      
      idx++;
      it++;
   }
   
   ii_Primes[idx] = 0;

   if (CpuSupportsAvx())
   {
      id_Primes[idx] = 0.0;
      SetMiniChunkRange(ip_PrimorialApp->GetMaxN() + 1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
   }
   
   ib_Initialized = true;
}

void  PrimorialWorker::CleanUp(void)
{
   xfree(ii_Primes);
   
   if (CpuSupportsAvx())
      xfree(id_Primes);
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
   
      if (primorial(ii_Primes, ps))
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
   // in this function.
   
   // We know that a factor was found, but the call to PrimorialWorker
   // doesn't tell us the term or the p for it, so we will use
   // the long way to find it.
   double      inverse, qd;
   uint32_t    maxN = ip_PrimorialApp->GetMaxN();
   uint32_t    value;
   uint32_t    term, idx;
   int64_t     q, r;
   long double pmax, vd;

   inverse = 1.0/p;
   r = 1;

   term = ii_Primes[0];
   idx = 0;
   value = 1;
   
   pmax = (long double) (p + 1);
   vd = (long double) value;   
   
   while (term > 0 && term <= maxN)
   {
      vd *= (long double) term;
      value *= term;

      // If the current value is greater than the
      // the max, then the value cannot be prime.      
      if (vd > pmax)
         break;
      
      if (value == p+1) 
         ip_PrimorialApp->ReportPrime(p, term, -1);
         
      if (value == p-1) 
         ip_PrimorialApp->ReportPrime(p, term, +1);
            
      idx++;
      term = ii_Primes[idx];
   }
   
   term = ii_Primes[0];
   idx = 0;
   
   while (term > 0 && term <= maxN)
   {
      qd = ((double) r * (double) term);

      q = (int64_t) (qd * inverse);

      r = (r*term) - p*q;

      if (r < 0)
         r += p;
      else if (r >= (int64_t) p)
         r -= p;
      
      if (r == +1) 
         ip_PrimorialApp->ReportFactor(p, term, -1);
      
      if (r == (int64_t) p-1) 
         ip_PrimorialApp->ReportFactor(p, term, +1);
      
      idx++;
      term = ii_Primes[idx];
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
   
   nextPrime[0] = 2.0;
   
   avx_set_1a(nextPrime);

   idx = 1;
   while (ii_Primes[idx] != 0)
   {
      nextPrime[0] = id_Primes[idx];
      
      avx_set_1b(nextPrime);
      avx_mulmod(dps, reciprocals);
      
      CheckAVXResult(miniPrimeChunk, dps, ii_Primes[idx]);
      
      idx++;
   }
}

void  PrimorialWorker::CheckAVXResult(uint64_t *ps, double *dps, uint32_t theN)
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

void  PrimorialWorker::VerifyAVXFactor(uint64_t p, uint32_t theN, int32_t theC)
{
   // Yes, I know that this doesn't check for primes like ExtractFactors, but we know
   // all primorial primes < 2^128, so it isn't important that this doesn't catch them.
   uint64_t rem = 1;
   uint32_t i = 0;
   bool     found = false;
   
   fpu_push_1divp(p);
   
   while (ii_Primes[i] <= theN)
   {     
      rem = fpu_mulmod(rem, ii_Primes[i], p);

      if (ii_Primes[i] == theN)
         break;
      
      i++;
   }
  
   fpu_pop();

   if (rem == +1)
   {
      ip_PrimorialApp->ReportFactor(p, theN, -1);
   
      if (theC == -1)
         found = true;
   }
   
   if (p == rem + 1)
   {
      ip_PrimorialApp->ReportFactor(p, theN, +1);
   
      if (theC == +1)
         found = true;
   }

 
   if (!found)
      FatalError("The AVX routine found factor %" PRIu64" of %u#%+d, but it could not be validated", p, theN, theC);
}