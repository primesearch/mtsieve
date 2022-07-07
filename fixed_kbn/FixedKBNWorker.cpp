/* FixedKBNWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "FixedKBNWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

FixedKBNWorker::FixedKBNWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_FixedKBNApp = (FixedKBNApp *) theApp;
   
   il_K = ip_FixedKBNApp->GetK();
   ii_Base = ip_FixedKBNApp->GetBase();
   ii_N = ip_FixedKBNApp->GetN();
   il_MinC = ip_FixedKBNApp->GetMinC();
   il_MaxC = ip_FixedKBNApp->GetMaxC();
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  FixedKBNWorker::CleanUp(void)
{
}

void  FixedKBNWorker::TestMegaPrimeChunk(void)
{
   uint64_t ps[4];
   uint64_t bs[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();   
      
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];

      bs[0] = bs[1] = bs[2] = bs[3] = ii_Base;
      
      fpu_powmod_4b_1n_4p(bs, ii_N, ps);
   
      RemoveTerms(ps[0], bs[0]);
      RemoveTerms(ps[1], bs[1]);
      RemoveTerms(ps[2], bs[2]);
      RemoveTerms(ps[3], bs[3]);

      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] > maxPrime)
         break;

      // If no terms left, then we are done
      if (ip_FixedKBNApp->GetTermCount() == 0) {
         SetLargestPrimeTested(maxPrime, 0);
         break;
      }
   }
}

void  FixedKBNWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("FixedKBNWorker::TestMiniPrimeChunk not implemented");
}

void    FixedKBNWorker::RemoveTerms(uint64_t prime, uint64_t bExpN)
{
   int64_t  c;
   int64_t  kbExpN;
   
   fpu_push_1divp(prime);
   
   kbExpN = fpu_mulmod(il_K, bExpN, prime);
   
   fpu_pop();
      
   // k*b^n (mod p) = kbExpN
   // k*b^n - kbExpN (mod p) = 0
   // k*b^n + (p - kbExpN) (mod p) = 0

   // Set c = p - kbExpN
   // Since p > kbExpN, c is positive
   // k*b^n + c (mod p) = 0
   c = (int64_t) (prime - kbExpN);

   int countVerified = 0;

   do 
   {
      // k*b^n + c (mod p) = 0
      ip_FixedKBNApp->ReportFactor(prime, c, (++countVerified < 5));

      c += (int64_t) prime;
   } while (c <= il_MaxC);
   
   c = (int64_t) (prime - kbExpN);

   countVerified = 0;
   do 
   {
      // k*b^n + c (mod p) = 0
      ip_FixedKBNApp->ReportFactor(prime, c, (++countVerified < 5));
      
      c -= (int64_t) prime;
   } while (c >= il_MinC);
}
