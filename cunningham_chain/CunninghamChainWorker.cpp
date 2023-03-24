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
   
   ib_HalfK = ip_CunninghamChainApp->IsHalfK();
   il_MinK = ip_CunninghamChainApp->GetMinK();
   il_MaxK = ip_CunninghamChainApp->GetMaxK();
   ii_Base = ip_CunninghamChainApp->GetBase();
   ii_N = ip_CunninghamChainApp->GetN();
   
   il_Terms = ip_CunninghamChainApp->GetTerms();

   ii_BaseInverses = NULL;
   il_MyPrimeList = NULL;
   ii_InverseList = NULL;
   
   // This limit comes from newpgen, but is not documented why it exists.
   // I do know that larger bases yield incorrect factors if we use the wrong logic.
   if (it_TermType == TT_BN && ii_Base > 2 && ii_Base < 255256)
      BuildBaseInverses();

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CunninghamChainWorker::CleanUp(void)
{
   if (il_MyPrimeList != NULL)
   {
      xfree(il_MyPrimeList);
      xfree(ii_InverseList);
   }
   
   if (ii_BaseInverses != NULL)
      xfree(ii_BaseInverses);
}

void  CunninghamChainWorker::NotifyPrimeListAllocated(uint32_t primesInList)
{
   if (il_MyPrimeList != NULL)
   {
      xfree(il_MyPrimeList);
      xfree(ii_InverseList);
   }
   
   if (it_TermType == TT_BN && ii_Base > 2 && ii_Base < 255256)
   {
      il_MyPrimeList = (uint64_t *) xmalloc((primesInList + 10) * sizeof(uint64_t));
      ii_InverseList = (uint32_t *) xmalloc((primesInList + 10) * sizeof(uint32_t));
   }
}

void  CunninghamChainWorker::TestMegaPrimeChunk(void)
{
   uint64_t ks[4];
   uint64_t ps[4];
   uint32_t idx;
   
   if (it_TermType == TT_BN && ii_Base > 2 && ii_Base < 255256)
   {
      TestSmallB();
      return;
   }

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
   }
}

void  CunninghamChainWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CunninghamChainWorker::TestMiniPrimeChunk not implemented");
}

void  CunninghamChainWorker::TestSmallB(void)
{
   uint64_t ks[4];
   uint64_t p1, ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t idx, count, svb, pmb;

   // Evaluate primes in the vector to determine if can yield a factor.  Only
   // put primes that can yield a factor into an array for the second loop.
   count = 0;
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx++)
   {
      p1 = il_PrimeList[pIdx];
         
      pmb = (p1 % ii_Base);
      
      if (ii_BaseInverses[pmb] == 0)
         continue;

      if (p1 > maxPrime)
         break;
      
      svb = ii_BaseInverses[pmb];
      
      il_MyPrimeList[count] = p1;
      ii_InverseList[count] = svb;
      count++;
   }
   
   if (count == 0)
      return;
   
   // Duplicate the last few entries so that the
   // number of valid entries is divisible by 4.
   while (count % 4 != 0)
   {
      il_MyPrimeList[count] = p1;
      ii_InverseList[count] = svb;
      count++;
   }

   for (idx=0; idx<count; idx+=4)
   {
      ps[0] = il_MyPrimeList[idx+0];
      ps[1] = il_MyPrimeList[idx+1];
      ps[2] = il_MyPrimeList[idx+2];
      ps[3] = il_MyPrimeList[idx+3];
      
      ks[0] = (1+ii_InverseList[idx+0]*ps[0])/ii_Base;
      ks[1] = (1+ii_InverseList[idx+1]*ps[1])/ii_Base;
      ks[2] = (1+ii_InverseList[idx+2]*ps[2])/ii_Base;
      ks[3] = (1+ii_InverseList[idx+3]*ps[3])/ii_Base;
      
      // Starting with k*2^n = 1 (mod p) 
      //           --> k = (1/2)^n (mod p)
      //           --> k = inverse^n (mod p)
      fpu_powmod_4b_1n_4p(ks, ii_N, &il_MyPrimeList[idx+0]);

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
   }

   // Adjust for the possibility that we tested the same prime
   // more than once at the end of the list.
   if (ps[3] == ps[2])
      SetLargestPrimeTested(ps[3], -1);
   if (ps[3] == ps[1])
      SetLargestPrimeTested(ps[3], -1);
   if (ps[3] == ps[0])
      SetLargestPrimeTested(ps[3], -1);
}

void  CunninghamChainWorker::RemoveTerms(uint64_t thePrime, uint64_t k, uint32_t termInChain)
{
   if (it_TermType == TT_BN)
   {
      // If b^n (mod p) == 0, then there is no solution for k.
      if (ii_Base >= thePrime && ii_Base % thePrime == 0)
         return;
   }
   else
   {
      // This can happen if n! (mod p) == 1 or n# (mod p) == 1 or n! (mod p) == -1 or n# (mod p) == -1
      if (thePrime == k || k == 0)
         return;
   }

   // Although k < thePrime, we need k >= il_MinK
   if (k < il_MinK)
   {
      if (thePrime >= il_MinK)
      {
         // We can avoid the use of division here
         k += thePrime; 
      }
      else
      {
         // Compute k such that il_MinK <= k < il_MinK + p
         k += thePrime * ((il_MinK - k + thePrime - 1)/thePrime);
      }
   }

   if (ib_HalfK)
   {
      if (ii_Base == 2)
      {
         // Ensure that k is odd
         if (!(k & 1))
            k += thePrime;
      }
      else
      {
         // Ensure that k is event
         if (k & 1)
            k += thePrime;
      }
   }

   if (k <= il_MaxK)
      ip_CunninghamChainApp->ReportFactor(thePrime, k, termInChain);
}

void    CunninghamChainWorker::BuildBaseInverses(void)
{
	ii_BaseInverses = (uint32_t *) xmalloc(ii_Base * sizeof(uint32_t));

	ii_BaseInverses[0] = 0;
	ii_BaseInverses[1] = ii_Base-1;

	for (uint32_t i = 2; i < ii_Base; i++)
		ii_BaseInverses[i] = EuclidExtendedGCD(i, ii_Base);
}

uint32_t CunninghamChainWorker::EuclidExtendedGCD(uint32_t a, uint32_t base)
{
	int64_t u = 0, d = base, v1 = 1, v3 = a;
   int64_t q, t1, t3;

	while (v3)
	{
		q = d/v3;
		t3 = d - (q*v3);
		t1 = u - (q*v1);

		u = v1;
		d = v3;
      
		v1 = t1;
		v3 = t3;
	}
   
	return (d == 1) ? (uint32_t)((u>0)?base-u:-u) : 0;
}