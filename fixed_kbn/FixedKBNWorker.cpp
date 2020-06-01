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
   uint64_t p1, p2, p3, p4, ps[4];
   uint64_t bs[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   
   vector<uint64_t>::iterator it = iv_Primes.begin(); 
   
   while (it != iv_Primes.end())
   {
      p1 = ps[0] = *it;
      it++;
      
      p2 = ps[1] = *it;
      it++;
      
      p3 = ps[2] = *it;
      it++;
      
      p4 = ps[3] = *it;
      it++;
      
      bs[0] = bs[1] = bs[2] = bs[3] = ii_Base;
      
      fpu_powmod_4b_1n_4p(bs, ii_N, ps);
   
      RemoveTerms(p1, bs[0]);
      RemoveTerms(p2, bs[1]);
      RemoveTerms(p3, bs[2]);
      RemoveTerms(p4, bs[3]);

      SetLargestPrimeTested(p4, 4);
      
      if (p4 > maxPrime)
         break;
   }
}

void  FixedKBNWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("FixedKBNWorker::TestMiniPrimeChunk not implemented");
}

void    FixedKBNWorker::RemoveTerms(uint64_t prime, uint64_t b)
{
   int64_t  c, cc;
   int64_t  sPrime;
   uint64_t rem;
   
   fpu_push_1divp(prime);
   
   c = prime - fpu_mulmod(il_K, b, prime);
   
   // Need this to be signed for the next calculation
   sPrime = (int64_t) prime;
   
   // Compute c such that il_MinC <= c < il_MinC + p
   cc = sPrime * ((il_MinC - c)/sPrime);

   c += cc;
   
   while (c > il_MaxC) c -= prime;
   while (c < il_MinC) c += prime;
   
   if (c > il_MaxC)
   {
      fpu_pop();
      return;
   }
      
   rem = fpu_powmod(ii_Base, ii_N, prime);
   il_KBpowN = fpu_mulmod(rem, il_K, prime);

   fpu_pop();
  
   do 
   {
      if (ip_FixedKBNApp->ReportFactor(prime, c))
         VerifyFactor(prime, c);

      c += prime;
   } while (c <= il_MaxC);
}


void  FixedKBNWorker::VerifyFactor(uint64_t prime, int64_t c)
{
   int64_t  rem;
   int64_t  sPrime;

   rem = c + il_KBpowN;
   
   // Need this to be signed for the next calculation
   sPrime = (int64_t) prime;
      
   rem = (rem % sPrime);

   if (rem < 0)
      rem += sPrime;
   
   if (rem != 0)
      FatalError("%" PRIu64"*%u^%u%+d mod %" PRIu64" = %" PRIu64"", il_K, ii_Base, ii_N, c, prime, rem);
}
