/* KBBWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "KBBWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"

KBBWorker::KBBWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_KBBApp = (KBBApp *) theApp;
   
   il_K = ip_KBBApp->GetK();

   // Give a little extra space
   ii_BaseCount = ip_KBBApp->GetMaxB() - ip_KBBApp->GetMinB() + 10;
   
   ip_Bases = (uint32_t *) xmalloc(ii_BaseCount * sizeof(uint32_t));
   
   il_NextBaseBuild = 0;
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  KBBWorker::CleanUp(void)
{
}

void  KBBWorker::TestMegaPrimeChunk(void)
{
   uint64_t p1, p2, p3, p4, ps[4];
   uint64_t bs[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t idx;
   
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
      
      // Every once in a while rebuild the list of base as it will have fewer entries
      // which will speed up testing for the next range of p.
      if (p1 > il_NextBaseBuild)
      {
         memset(ip_Bases, 0, ii_BaseCount * sizeof(uint32_t));
         
         ip_KBBApp->GetBases(ip_Bases);
         
         il_NextBaseBuild = (p4 << 1);
      }
   
      for (idx=0; idx<ii_BaseCount; idx++)
      {
         if (ip_Bases[idx] == 0)
            break;
         
         bs[0] = bs[1] = bs[2] = bs[3] = ip_Bases[idx];
         
         sse_powmod_4b_1n_4p_mulmod_1k(bs,  bs[0], ps, il_K);
         
         // bs[0] = il_K*b^b mod p1
         CheckK(p1, ip_Bases[idx], bs[0]);
         
         // bs[1] = il_K*b^b mod p2
         CheckK(p2, ip_Bases[idx], bs[1]);
         
         // bs[2] = il_K*b^b mod p3
         CheckK(p3, ip_Bases[idx], bs[2]);
         
         // bs[3] = il_K*b^b mod p3
         CheckK(p4, ip_Bases[idx], bs[3]);
      }

      SetLargestPrimeTested(p4, 4);
      
      if (p4 > maxPrime)
         break;
   }
}

void  KBBWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("KBBWorker::TestMiniPrimeChunk not implemented");
}

void  KBBWorker::CheckK(uint64_t prime, uint64_t base, uint64_t rem)
{   
   if (rem == +1)
   {
      if (ip_KBBApp->ReportFactor(prime, base, -1))
         VerifyFactor(prime, base, -1);
   }

   if (rem == prime-1)
   {
      if (ip_KBBApp->ReportFactor(prime, base, +1))
         VerifyFactor(prime, base, +1);
   }
}
   
void  KBBWorker::VerifyFactor(uint64_t prime, uint64_t b, int32_t c)
{
   uint64_t  rem;
   char      bStr[50];
   char      kStr[50];
   char      pStr[50];
   char      rStr[50];

   fpu_push_1divp(prime);
   
   rem = fpu_powmod(b, b, prime);
   rem = fpu_mulmod(rem, il_K, prime);

   fpu_pop();
   
   if (c == +1 && rem != prime-1)
   {
      sprintf(bStr, "%" PRIu64"", b);
      sprintf(kStr, "%" PRIu64"", il_K);
      sprintf(pStr, "%" PRIu64"", prime);
      sprintf(rStr, "%" PRIu64"", rem);
      
      FatalError("%s*%s^%s+1 mod %s = %s", kStr, bStr, bStr, pStr, rStr);
   }
   
   if (c == -1 && rem != +1)
   {
      sprintf(bStr, "%" PRIu64"", b);
      sprintf(kStr, "%" PRIu64"", il_K);
      sprintf(pStr, "%" PRIu64"", prime);
      sprintf(rStr, "%" PRIu64"", rem);
      
      FatalError("%s*%s^%s-1 mod %s = %s", kStr, bStr, bStr, pStr, rStr);
   }
}
