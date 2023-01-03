/* CunninghamChainWorker.cpp -- (C) Mark Rodenkirch, June 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "../core/MpArithVector.h"
#include "../x86_asm/fpu-asm-x86.h"

#include "CunninghamChainWorker.h"

CunninghamChainWorker::CunninghamChainWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_CunninghamChainApp = (CunninghamChainApp *) theApp;
   
   it_ChainKind = ip_CunninghamChainApp->GetChainKind();
   it_TermType = ip_CunninghamChainApp->GetTermType();
   ii_ChainLength = ip_CunninghamChainApp->GetChainLength();
   
   il_MinK = ip_CunninghamChainApp->GetMinK();
   il_MaxK = ip_CunninghamChainApp->GetMaxK();
   ii_Base = ip_CunninghamChainApp->GetBase();
   ii_N = ip_CunninghamChainApp->GetN();
   
   il_Terms = ip_CunninghamChainApp->GetTerms();

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CunninghamChainWorker::CleanUp(void)
{
}

void  CunninghamChainWorker::TestMegaPrimeChunk(void)
{
   uint64_t ks[4];
   uint64_t ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t idx;
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      if (it_TermType == TT_BN)
      {
         // Compute ks as (1/b) (mod p)
         if (ii_Base == 2)
         {
            ks[0] = (1+ps[0]) >> 1;
            ks[1] = (1+ps[1]) >> 1;
            ks[2] = (1+ps[2]) >> 1;
            ks[3] = (1+ps[3]) >> 1;
         }
         else
         {
            ks[0] = InvMod32(ii_Base, ps[0]);
            ks[1] = InvMod32(ii_Base, ps[1]);
            ks[2] = InvMod32(ii_Base, ps[2]);
            ks[3] = InvMod32(ii_Base, ps[3]);
         }
         
         // ks = (1/b)^n (mod p)
         fpu_powmod_4b_1n_4p(ks, ii_N, ps);
      }
      else
      {
         MpArithVec mp(ps);

         MpResVec resRem = mp.one();
      
         idx = 0;
         while (il_Terms[idx] > 0)
         {
            resRem = mp.mul(resRem, mp.nToRes(il_Terms[idx]));

            idx++;
         }
      
         resRem = mp.resToN(resRem);

         ks[0] = InvMod64(resRem[0], ps[0]);
         ks[1] = InvMod64(resRem[1], ps[1]);
         ks[2] = InvMod64(resRem[2], ps[2]);
         ks[3] = InvMod64(resRem[3], ps[3]);
      }
      
      if (it_ChainKind == CCT_SECONDKIND)
      {
         ks[0] = ps[0] - ks[0];
         ks[1] = ps[1] - ks[1];
         ks[2] = ps[2] - ks[2];
         ks[3] = ps[3] - ks[3];
      }
      
      // In the first iteration of this loop we have computed k such that k*m# (mod p) = +1 or -1.
      // For subsequent iterations we are computing k such that 2^n*k*m# (mod p) = +1 or -1.
      // In short p divides the first term of the Cunningham Chain for k1 and
      // p divides the second term of the Cunningham Chain for k2. 
      for (uint32_t termInChain=1; termInChain<=ii_ChainLength; termInChain++)
      {
         // Now we have k such that:
         //    k*m#+1 (mod p) = 0 (first kind)
         // or k*m#-1 (mod p) = 0 (second kind)

         // Note that 0 <= k < p
         RemoveTerms(ps[0], ks[0], termInChain);
         RemoveTerms(ps[1], ks[1], termInChain);
         RemoveTerms(ps[2], ks[2], termInChain);
         RemoveTerms(ps[3], ks[3], termInChain);

         // Make sure k is even before dividing by 2.
         if (ks[0] & 1) ks[0] += ps[0];
         if (ks[1] & 1) ks[1] += ps[1];
         if (ks[2] & 1) ks[2] += ps[2];
         if (ks[3] & 1) ks[3] += ps[3];
            
         ks[0] >>= 1;
         ks[1] >>= 1;
         ks[2] >>= 1;
         ks[3] >>= 1;
      }

      SetLargestPrimeTested(ps[3], 4);
   
      if (ps[3] >= maxPrime)
         break;
   }
}

void  CunninghamChainWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CunninghamChainWorker::TestMiniPrimeChunk not implemented");
}

void    CunninghamChainWorker::RemoveTerms(uint64_t prime, uint64_t k, uint32_t termInChain)
{
   if (it_TermType == TT_BN)
   {
      // If b^n (mod p) == 0, then there is no solution for k.
      if (ii_Base >= prime && ii_Base % prime == 0)
         return;
   }
   else
   {
      // This can happen if n! (mod p) == 1 or n# (mod p) == 1 or n! (mod p) == -1 or n# (mod p) == -1
      if (prime == k || k == 0)
         return;
   }

   // Although k < prime, we need k >= il_MinK
   if (k < il_MinK)
   {
      if (prime >= il_MinK)
      {
         // We can avoid the use of division here
         k += prime; 
      }
      else
      {
         // Compute k such that il_MinK <= k < il_MinK + p
         k += prime * ((il_MinK - k + prime - 1)/prime);
      }
   }

   if (k <= il_MaxK)
      ip_CunninghamChainApp->ReportFactor(prime, k, termInChain);
}
