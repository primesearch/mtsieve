/* SmarandacheWellinWorker.cpp -- (C) Mark Rodenkirch, March 2025

   This program is free software;you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation;either version 2 of the License, or
   (at your option) any later version.
*/

#include <math.h>
#include "SmarandacheWellinWorker.h"
#include "../core/MpArithVector.h"

extern "C" int mfsieve(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);
extern "C" int SmarandacheWellin(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);

SmarandacheWellinWorker::SmarandacheWellinWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_SmarandacheWellinApp = (SmarandacheWellinApp *) theApp;
   
   ii_MinN = ip_SmarandacheWellinApp->GetMinN();
   ii_MaxN = ip_SmarandacheWellinApp->GetMaxN();

   ip_Primes = ip_SmarandacheWellinApp->GetPrimes(ii_NumberOfPrimes);
   
   ib_Initialized = true;
}

void  SmarandacheWellinWorker::CleanUp(void)
{
}

void  SmarandacheWellinWorker::TestMegaPrimeChunk(void)
{
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   uint32_t  i;
     
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];

      MpArithVec mp(ps);
      MpResVec   zero = mp.zero();
      MpResVec   res = mp.nToRes(2357);
      MpResVec   mp1e2 = mp.nToRes(100);
      MpResVec   mp1e3 = mp.nToRes(1000);
      MpResVec   mp1e4 = mp.nToRes(10000);
      MpResVec   mp1e5 = mp.nToRes(100000);
      MpResVec   mp1e6 = mp.nToRes(1000000);
      MpResVec   mp1e7 = mp.nToRes(10000000);
      MpResVec   mp1e8 = mp.nToRes(100000000);
      MpResVec   mp1e9 = mp.nToRes(1000000000);

      for (i=0; i<ii_NumberOfPrimes; i++)
      {
         if (ip_Primes[i] < 100)
            res = mp.mul(res, mp1e2);
         else if (ip_Primes[i] < 1000)
            res = mp.mul(res, mp1e3);
         else if (ip_Primes[i] < 10000)
            res = mp.mul(res, mp1e4);
         else if (ip_Primes[i] < 100000)
            res = mp.mul(res, mp1e5);
         else if (ip_Primes[i] < 1000000)
            res = mp.mul(res, mp1e6);
         else if (ip_Primes[i] < 10000000)
            res = mp.mul(res, mp1e7);
         else if (ip_Primes[i] < 100000000)
            res = mp.mul(res, mp1e8);
         else if (ip_Primes[i] < 1000000000)
            res = mp.mul(res, mp1e9);
         else
            FatalError("How did this happen?  n is too large");
         
         res = mp.add(res, mp.nToRes(ip_Primes[i]));

         if (res[0] == zero[0])
            ip_SmarandacheWellinApp->ReportFactor(ps[0], ip_Primes[i]);
         if (res[1] == zero[1])
            ip_SmarandacheWellinApp->ReportFactor(ps[1], ip_Primes[i]);
         if (res[2] == zero[2])
            ip_SmarandacheWellinApp->ReportFactor(ps[2], ip_Primes[i]);
         if (res[3] == zero[3])
            ip_SmarandacheWellinApp->ReportFactor(ps[3], ip_Primes[i]);         
      }
      
      SetLargestPrimeTested(ps[3], 4);
   
      if (ps[3] >= maxPrime)
         break;
   }
}

void  SmarandacheWellinWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheWellinWorker::TestMiniPrimeChunk not implemented");
}
