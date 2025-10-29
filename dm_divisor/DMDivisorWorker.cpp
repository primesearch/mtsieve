/* DMDivisorWorker.cpp -- (C) Mark Rodenkirch, September 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "DMDivisorWorker.h"
#include "../core/MpArithVector.h"

DMDivisorWorker::DMDivisorWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_DMDivisorApp = (DMDivisorApp *) theApp;
   
   il_MinK = ip_DMDivisorApp->GetMinK();
   il_MaxK = ip_DMDivisorApp->GetMaxK();
   ii_N = ip_DMDivisorApp->GetN();
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  DMDivisorWorker::CleanUp(void)
{
}

void  DMDivisorWorker::TestMegaPrimeChunk(void)
{
   uint64_t k1, k2, k3, k4;
   uint64_t bs[4], ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
            
      bs[0] = bs[1] = bs[2] = bs[3] = 2;
      
      MpArithVec mp(ps);
      MpResVec   two = mp.nToRes(2);
      MpResVec   res = mp.pow(two, ii_N);
      
      res = mp.resToN(res);

      // Now bs = 2^exp-1 (mod p)
      bs[0] = res[0] - 1;
      bs[1] = res[1] - 1;
      bs[2] = res[2] - 1;
      bs[3] = res[3] - 1;
      
      // We are looking for k such that 2*k*bs+1 (mod p) = 0
      k1 = InvMod64(bs[0], ps[0]);
      k2 = InvMod64(bs[1], ps[1]);
      k3 = InvMod64(bs[2], ps[2]);
      k4 = InvMod64(bs[3], ps[3]);
      
      // 2*k*bs+1 = 0 (mod p) --> 2*k = -invbs (mod p)
      
      // Now ensure that invbs is positive
      k1 = ps[0] - k1;
      k2 = ps[1] - k2;
      k3 = ps[2] - k3;
      k4 = ps[3] - k4;
      
      // We need invbs to be even so that we can divide by 2
      if (k1 & 1) k1 += ps[0];
      if (k2 & 1) k2 += ps[1];
      if (k3 & 1) k3 += ps[2];
      if (k4 & 1) k4 += ps[3];
      
      k1 >>= 1;
      k2 >>= 1;
      k3 >>= 1;
      k4 >>= 1;
      
      // We have now solved for mink
      if (k1 <= il_MaxK) RemoveTerms(ps[0], k1);
      if (k2 <= il_MaxK) RemoveTerms(ps[1], k2);
      if (k3 <= il_MaxK) RemoveTerms(ps[2], k3);
      if (k4 <= il_MaxK) RemoveTerms(ps[3], k4);

      SetLargestPrimeTested(ps[3], 4);
   
      if (ps[3] >= maxPrime)
         break;

      if (ip_App->IsInterrupted())
         break;
   }
}

void  DMDivisorWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("DMDivisorWorker::TestMiniPrimeChunk not implemented");
}

// This must be used when prime <= il_MaxK.
void    DMDivisorWorker::RemoveTerms(uint64_t prime, uint64_t k)
{
   // Adjust so that k >= il_MinK

   if (k < il_MinK)
   {
      if (prime >= il_MinK)
         k += prime; 
      else
      {
         // Compute k such that il_MinK <= k < il_MinK + p
         k += prime * ((il_MinK - k + prime - 1)/prime);
      }
   }
   
   if (k > il_MaxK)
      return;
   
   ip_DMDivisorApp->ReportFactor(prime, k);
}
