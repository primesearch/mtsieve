/* AlternatingFactorial.cpp -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "AlternatingFactorialWorker.h"

extern "C" {
   int afsieve(const uint32_t nmax, const uint64_t *P);
}

AlternatingFactorialWorker::AlternatingFactorialWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_AlternatingFactorialApp = (AlternatingFactorialApp *) theApp;
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  AlternatingFactorialWorker::CleanUp(void)
{
}

void  AlternatingFactorialWorker::TestMegaPrimeChunk(void)
{
   uint64_t  p[4], maxPrime = ip_App->GetMaxPrime();
   uint32_t  gotFactor, maxN = ip_AlternatingFactorialApp->GetMaxN();
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
      
      gotFactor = afsieve(maxN, p);

      if (gotFactor)
      {
         ExtractFactors(p[0]);
         ExtractFactors(p[1]);
         ExtractFactors(p[2]);
         ExtractFactors(p[3]);
      }
      
      SetLargestPrimeTested(p[3], 4);
      
      if (p[3] >= maxPrime)
         break;
   }
}

void  AlternatingFactorialWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("AlternatingFactorialWorker::TestMiniPrimeChunk not implemented");
}

void  AlternatingFactorialWorker::ExtractFactors(uint64_t p)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.

   // We know that a factor was found, but the call to AlternatingFactorial
   // doesn't tell us the term or the p for it, so we will use
   // the long way to find it.
   double   inverse, qd;
   uint32_t n;

   int64_t  up = p;
   int64_t  q, rem;
   int64_t  nm1term;
   
   inverse = 1.0 / p;
   nm1term = 1;
   rem = 1;
   
   for (n=2; n<=ip_AlternatingFactorialApp->GetMaxN(); n++)
   {
      qd = ((double) rem * (double) n);

      q = (int64_t) (qd * inverse);

      rem = (rem*n) - p*q;

      if (rem < 0)
         rem += p;
      else if (rem >= up)
         rem -= p;
      
      // rem = n! % p
      // nm1term = af(n-1) % p
      if (rem == nm1term)
         ip_AlternatingFactorialApp->ReportFactor(p, n);

      // Set nm1term to af(n) % p
      if (rem > nm1term)
         nm1term = rem - nm1term;
      else
         nm1term = (rem + p - nm1term);
      
   }
}
