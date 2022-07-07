/* SophieGermainWorker.cpp -- (C) Mark Rodenkirch, July 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "../core/MpArithVector.h"

#include "SophieGermainWorker.h"

SophieGermainWorker::SophieGermainWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_SophieGermainApp = (SophieGermainApp *) theApp;
   
   il_MinK = ip_SophieGermainApp->GetMinK();
   il_MaxK = ip_SophieGermainApp->GetMaxK();
   ii_Base = ip_SophieGermainApp->GetBase();
   ii_N = ip_SophieGermainApp->GetN();
   ib_GeneralizedSearch = ip_SophieGermainApp->IsGeneralizedSearch();
     
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  SophieGermainWorker::CleanUp(void)
{
}

void  SophieGermainWorker::TestMegaPrimeChunk(void)
{
   if (il_PrimeList[0] < il_MaxK || il_PrimeList[0] < 50)
      TestMegaPrimeChunkSmall();
   else
      TestMegaPrimeChunkLarge();
}

void  SophieGermainWorker::TestMegaPrimeChunkSmall(void)
{
   uint64_t invs[4];
   uint64_t ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();

   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];

      if (ii_Base == 2)
      {
         invs[0] = (1+ps[0]) >> 1;
         invs[1] = (1+ps[1]) >> 1;
         invs[2] = (1+ps[2]) >> 1;
         invs[3] = (1+ps[3]) >> 1;
      }
      else
      {
         invs[0] = InvMod32(ii_Base, ps[0]);
         invs[1] = InvMod32(ii_Base, ps[1]);
         invs[2] = InvMod32(ii_Base, ps[2]);
         invs[3] = InvMod32(ii_Base, ps[3]);
      }
      
      // Starting with k*b^n = 1 (mod p) 
      //           --> k = (1/b)^n (mod p)
      //           --> k = inv(b)^n (mod p)
      MpArithVec mp(ps);
      
      MpResVec resInvs = mp.nToRes(invs);
      MpResVec res = mp.pow(resInvs, ii_N);
      MpResVec resKs = mp.resToN(res);

      if (resKs[0] <= il_MaxK) RemoveTermsSmallPrime(resKs[0], true, ps[0]);
      if (resKs[1] <= il_MaxK) RemoveTermsSmallPrime(resKs[1], true, ps[1]);
      if (resKs[2] <= il_MaxK) RemoveTermsSmallPrime(resKs[2], true, ps[2]);
      if (resKs[3] <= il_MaxK) RemoveTermsSmallPrime(resKs[3], true, ps[3]);

      if (ii_Base == 2)
      {
         res = mp.mul(res, resInvs);
      }
      else
      {
         if (ib_GeneralizedSearch)
         {
            // Multipley by inv(ii_Base)
            res = mp.mul(res, resInvs);
         }
         else
         {
            // Multipley by inv(2)
            invs[0] = (1+ps[0]) >> 1;
            invs[1] = (1+ps[1]) >> 1;
            invs[2] = (1+ps[2]) >> 1;
            invs[3] = (1+ps[3]) >> 1;
            
            res = mp.mul(res, mp.nToRes(invs));
         }
      }
      
      resKs = mp.resToN(res);
      
      if (resKs[0] <= il_MaxK) RemoveTermsSmallPrime(resKs[0], false, ps[0]);
      if (resKs[1] <= il_MaxK) RemoveTermsSmallPrime(resKs[1], false, ps[1]);
      if (resKs[2] <= il_MaxK) RemoveTermsSmallPrime(resKs[2], false, ps[2]);
      if (resKs[3] <= il_MaxK) RemoveTermsSmallPrime(resKs[3], false, ps[3]);
      
      SetLargestPrimeTested(ps[3], 4);
   
      if (ps[3] >= maxPrime)
         break;
   }
}

void  SophieGermainWorker::TestMegaPrimeChunkLarge(void)
{
   uint64_t invs[4];
   uint64_t ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t pIdx = 0;

   while (pIdx < ii_PrimesInList)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      pIdx += 4;

      if (ii_Base == 2)
      {
         invs[0] = (1+ps[0]) >> 1;
         invs[1] = (1+ps[1]) >> 1;
         invs[2] = (1+ps[2]) >> 1;
         invs[3] = (1+ps[3]) >> 1;
      }
      else
      {
         invs[0] = InvMod32(ii_Base, ps[0]);
         invs[1] = InvMod32(ii_Base, ps[1]);
         invs[2] = InvMod32(ii_Base, ps[2]);
         invs[3] = InvMod32(ii_Base, ps[3]);
      }

      // Starting with k*b^n = 1 (mod p) 
      //           --> k = (1/b)^n (mod p)
      //           --> k = inv(b)^n (mod p)
      MpArithVec mp(ps);

      MpResVec resInvs = mp.nToRes(invs);
      MpResVec res = mp.pow(resInvs, ii_N);
      MpResVec resKs = mp.resToN(res);

      if (resKs[0] >= il_MinK && resKs[0] <= il_MaxK) RemoveTermsLargePrime(resKs[0], true, ps[0]);
      if (resKs[1] >= il_MinK && resKs[1] <= il_MaxK) RemoveTermsLargePrime(resKs[1], true, ps[1]);
      if (resKs[2] >= il_MinK && resKs[2] <= il_MaxK) RemoveTermsLargePrime(resKs[2], true, ps[2]);
      if (resKs[3] >= il_MinK && resKs[3] <= il_MaxK) RemoveTermsLargePrime(resKs[3], true, ps[3]);
      
      if (ii_Base == 2)
      {
         res = mp.mul(res, resInvs);
      }
      else 
      {
         invs[0] = (1+ps[0]) >> 1;
         invs[1] = (1+ps[1]) >> 1;
         invs[2] = (1+ps[2]) >> 1;
         invs[3] = (1+ps[3]) >> 1;
         res = mp.mul(res, mp.nToRes(invs));
      }
      
      resKs = mp.resToN(res);
      
      RemoveTermsLargePrime(resKs[0], false, ps[0]);
      RemoveTermsLargePrime(resKs[1], false, ps[1]);
      RemoveTermsLargePrime(resKs[2], false, ps[2]);
      RemoveTermsLargePrime(resKs[3], false, ps[3]);

      SetLargestPrimeTested(ps[3], 4);
   
      if (ps[3] >= maxPrime)
         break;
   }
}

void  SophieGermainWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SophieGermainWorker::TestMiniPrimeChunk not implemented");
}

// This must be used when prime <= il_MaxK.
void    SophieGermainWorker::RemoveTermsSmallPrime(uint64_t k, bool firstOfPair, uint64_t prime)
{
   // Make sure that k >= il_MinK
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

   if (ii_Base & 1)
   {
      // Make sure that k is even
      if (k & 1)
         k += prime;
   }
   else
   {
      // Make sure that k is odd
      if (!(k & 1))
         k += prime;
   }

   if (k > il_MaxK)
      return;

   uint32_t verifiedCount = 0;
   
   do
	{
      verifiedCount++;
      
      // If the first few are valid, then assume that the rest are valid.  This will
      // speed up testing of small primes.
      ip_SophieGermainApp->ReportFactor(prime, k, firstOfPair, (verifiedCount < 5));

		k += (prime << 1); 
	} while (k <= il_MaxK);
}


// Using this bypasses a number of if checks that can be done when prime > il_MaxK.
// Do not report k/n combinations if the k*2^n+1 is divisible by any p < 50
void    SophieGermainWorker::RemoveTermsLargePrime(uint64_t k, bool firstOfPair, uint64_t prime)
{   
   uint32_t smallN;
   uint64_t smallK;

   if (ii_Base == 2)
   {
      smallN = ii_N % (2);
      smallK = k % (3);
      if ((smallK << smallN) % (3) == 2) return;
        
      smallN = ii_N % (4);
      smallK = k % (5);
      if ((smallK << smallN) % (5) == 4) return;
      
      smallN = ii_N % (6);
      smallK = k % (7);
      if ((smallK << smallN) % (7) == 6) return;
      
      smallN = ii_N % (10);
      smallK = k % (11);
      if ((smallK << smallN) % (11) == 10) return;
      
      smallN = ii_N % (12);
      smallK = k % (13);
      if ((smallK << smallN) % (13) == 12) return;
      
      smallN = ii_N % (16);
      smallK = k % (17);
      if ((smallK << smallN) % (17) == 16) return;
      
      smallN = ii_N % (18);
      smallK = k % (19);
      if ((smallK << smallN) % (19) == 18) return;
      
      smallN = ii_N % (22);
      smallK = k % (23);
      if ((smallK << smallN) % (23) == 22) return;
      
      smallN = ii_N % (28);
      smallK = k % (29);
      if ((smallK << smallN) % (29) == 28) return;
      
      smallN = ii_N % (30);
      smallK = k % (31);
      if ((smallK << smallN) % (31) == 30) return;

      smallN = ii_N % (36);
      smallK = k % (37);
      if ((smallK << smallN) % (37) == 36) return;
      
      smallN = ii_N % (40);
      smallK = k % (41);
      if ((smallK << smallN) % (41) == 40) return;
      
      smallN = ii_N % (42);
      smallK = k % (43);
      if ((smallK << smallN) % (43) == 42) return;
      
      smallN = ii_N % (46);
      smallK = k % (47);
      if ((smallK << smallN) % (47) == 46) return;
   }

   if (ii_Base & 1)
   {
      // Make sure that k is even
      if (k & 1)
         k += prime;
   }
   else
   {
      // Make sure that k is odd
      if (!(k & 1))
         k += prime;
   }

   if (k < il_MinK || k > il_MaxK) return;

   ip_SophieGermainApp->ReportFactor(prime, k, firstOfPair, true);
}
