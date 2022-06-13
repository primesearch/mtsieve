/* PrimesInXWorker.cpp -- (C) Mark Rodenkirch, October 2012

   This class sets up the call to the PrimesInXWorker GPU function and parses the output
   from the GPU to determine if we have a PrimesInXWorker prime.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "PrimesInXWorker.h"

extern "C" int   pixsieve(const uint32_t *N1, const uint32_t *N2, const uint64_t *P, const uint64_t mult) __attribute__ ((pure));

PrimesInXWorker::PrimesInXWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_PrimesInXApp = (PrimesInXApp *) theApp;

   // Each thread needs its own copy of the terms
   ii_e3Terms = ip_PrimesInXApp->Get3DigitTermsCopy();
   ii_e6Terms = ip_PrimesInXApp->Get6DigitTermsCopy();
   ii_e9Terms = ip_PrimesInXApp->Get9DigitTermsCopy();
   
   ib_Initialized = true;
}

void  PrimesInXWorker::CleanUp(void)
{   
   if (ii_e3Terms) xfree(ii_e3Terms);
   if (ii_e6Terms) xfree(ii_e6Terms);
   if (ii_e9Terms) xfree(ii_e9Terms);
}

void  PrimesInXWorker::TestMegaPrimeChunk(void)
{
   uint64_t  ps[4];
   uint32_t  saveTerm;
   int32_t   gotFactor, index = 0, power;
   uint32_t  lminRemaining = ip_PrimesInXApp->GetMinLengthRemaining();
   uint32_t *n1, *n2 = ip_PrimesInXApp->Get1DigitTerms();
   uint32_t  multiplier;
   uint64_t  maxPrime = ip_PrimesInXApp->GetMaxPrime();

   for (uint32_t pIdx=0; pIdx<ii_WorkSize; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
            
      if (lminRemaining > 9 && ps[0] > 1000000000)
      {
         power = 9;
         
         // Can only use this after reaching p=1e9
         index = (lminRemaining / power) - 1;
         n1 = ii_e9Terms;
         multiplier = 1000000000;
      }
      else if (lminRemaining > 6 && ps[0] > 1000000)
      {
         power = 6;
         
         // Can only use this after reaching p=1e6
         index = (lminRemaining / power) - 1;
         n1 = ii_e6Terms;
         multiplier = 1000000;
      }
      else if (lminRemaining > 3 && ps[0] > 1000)
      {
         power = 3;
         
         // Can only use this after reaching p=1e3
         index = (lminRemaining / power) - 1;
         n1 = ii_e3Terms;
         multiplier = 1000;
      }

      if (index > 0)
      {
         // Set the "end of list" dynamically
         saveTerm = n1[index];
         n1[index] = multiplier;
      
         gotFactor = pixsieve(n1, &n2[index*power], ps, multiplier);
      
         // Restore the "end of list"
         n1[index] = saveTerm;
      }
      else
         gotFactor = pixsieve(n2, 0, ps, 10);
   
      if (gotFactor)
      {
         ExtractFactors(ps[0]);
         ExtractFactors(ps[1]);
         ExtractFactors(ps[2]);
         ExtractFactors(ps[3]);
      }
      
      SetLargestPrimeTested(ps[3], 4);

      if (ps[0] > maxPrime)
         break;
   }
}

void  PrimesInXWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("PrimesInXWorker::TestMiniPrimeChunk not implemented");
}

void  PrimesInXWorker::ExtractFactors(uint64_t p)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.

   // We know that a factor was found, but the call to pixsieve
   // doesn't tell us the term or the p for it, so we will use
   // the long way to find it.
   double   inverse, qd;
   uint64_t rem;
   uint32_t i, n, plength, firstDigit;
   int64_t  q;

   uint32_t *terms = ip_PrimesInXApp->Get1DigitTerms();

   plength = 1;
   rem = p;
   while (rem > 10)
   {
      plength++;
      rem /= 10;
   }
   firstDigit = rem;
   
   inverse = 1.0/p;
   rem = 0;

   for (i=0; i<=ip_PrimesInXApp->GetMaxLength(); i++)
   {
#if defined(__GNUC__) && defined(PREFETCH)
      __builtin_prefetch(&terms[i+PREFETCH], 0, 0);
#endif
      n = terms[i];
      
      qd = (((double)rem * 10.0) + (double)n);

      q = (int64_t) (qd * inverse);

      rem = (rem*10)+n - p*q;

      if (rem < 0)
         rem += p;
      else if (rem >= p)
         rem -= p;

      if (rem == 0)
      {
         if (plength == i+1 && terms[0] == firstDigit)
            ip_PrimesInXApp->ReportPrime(p, i+1);
         else 
            ip_PrimesInXApp->ReportFactor(p, i+1);
      }
   }
}
