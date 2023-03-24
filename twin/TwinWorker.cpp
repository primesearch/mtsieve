/* TwinWorker.cpp -- (C) Mark Rodenkirch, September 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "TwinWorker.h"
#include "../core/MpArithVector.h"
#include "../x86_asm/fpu-asm-x86.h"

TwinWorker::TwinWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_TwinApp = (TwinApp *) theApp;
   
   il_MinK = ip_TwinApp->GetMinK();
   il_MaxK = ip_TwinApp->GetMaxK();
   ii_Base = ip_TwinApp->GetBase();
   ii_N = ip_TwinApp->GetN();
   it_TermType = ip_TwinApp->GetTermType();
   
   ii_BaseInverses = NULL;
   il_MyPrimeList = NULL;
   ii_InverseList = NULL;
   
   il_Terms = ip_TwinApp->GetTerms();
   
   DeterminePrimeTermRange();
   
   // This limit comes from newpgen, but is not documented why it exists.
   // I do know that larger bases yield incorrect factors if we use the wrong logic.
   if (ii_Base < 255256 && it_TermType == TT_BN)
      BuildBaseInverses();
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  TwinWorker::CleanUp(void)
{
   if (il_MyPrimeList != NULL)
   {
      xfree(il_MyPrimeList);
      xfree(ii_InverseList);
   }
   
   if (ii_BaseInverses != NULL)
      xfree(ii_BaseInverses);
}

void  TwinWorker::NotifyPrimeListAllocated(uint32_t primesInList)
{
   if (il_MyPrimeList != NULL)
   {
      xfree(il_MyPrimeList);
      xfree(ii_InverseList);
   }
   
   if (ii_Base < 255256 && it_TermType == TT_BN)
   {
      il_MyPrimeList = (uint64_t *) xmalloc((primesInList + 10) * sizeof(uint64_t));
      ii_InverseList = (uint32_t *) xmalloc((primesInList + 10) * sizeof(uint32_t));
   }
}

void  TwinWorker::TestMegaPrimeChunk(void)
{
   if (ii_Base < 255256 && it_TermType == TT_BN)
      TestSmallB();
   else
      TestLargeB();
}


void  TwinWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("TwinWorker::TestMiniPrimeChunk not implemented");
}

void  TwinWorker::TestSmallB(void)
{
   uint64_t p1 = 0, p2, p3, p4;
   uint64_t k1, k2, k3, k4, ks[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   int32_t  svb = 0;
   int32_t  pmb, count, idx;
   
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
      p1 = il_MyPrimeList[idx+0];
      p2 = il_MyPrimeList[idx+1];
      p3 = il_MyPrimeList[idx+2];
      p4 = il_MyPrimeList[idx+3];
      
      ks[0] = k1 = (1+ii_InverseList[idx+0]*p1)/ii_Base;
      ks[1] = k2 = (1+ii_InverseList[idx+1]*p2)/ii_Base;
      ks[2] = k3 = (1+ii_InverseList[idx+2]*p3)/ii_Base;
      ks[3] = k4 = (1+ii_InverseList[idx+3]*p4)/ii_Base;
      
      // Starting with k*2^n = 1 (mod p) 
      //           --> k = (1/2)^n (mod p)
      //           --> k = inverse^n (mod p)
      fpu_powmod_4b_1n_4p(ks, ii_N, &il_MyPrimeList[idx+0]);
            
      k1 = ks[0];
      k2 = ks[1];
      k3 = ks[2];
      k4 = ks[3];

      if (p1 <= il_MaxK)
      {
         // This function allows for multiple terms to be removed
         if (k1 <= il_MaxK) RemoveTermsSmallPrime(p1, k1, -1);
         if (k2 <= il_MaxK) RemoveTermsSmallPrime(p2, k2, -1);
         if (k3 <= il_MaxK) RemoveTermsSmallPrime(p3, k3, -1);
         if (k4 <= il_MaxK) RemoveTermsSmallPrime(p4, k4, -1);
      }
      else
      {
         // This function allows for a single term to be removed
         if (k1 <= il_MaxK) RemoveTermsBigPrime(p1, k1, -1);
         if (k2 <= il_MaxK) RemoveTermsBigPrime(p2, k2, -1);
         if (k3 <= il_MaxK) RemoveTermsBigPrime(p3, k3, -1);
         if (k4 <= il_MaxK) RemoveTermsBigPrime(p4, k4, -1);
      }
      
      k1 = p1 - ks[0];
      k2 = p2 - ks[1];
      k3 = p3 - ks[2];
      k4 = p4 - ks[3];
         
      if (p1 <= il_MaxK)
      {
         // This function allows for multiple terms to be removed
         if (k1 <= il_MaxK) RemoveTermsSmallPrime(p1, k1, +1);
         if (k2 <= il_MaxK) RemoveTermsSmallPrime(p2, k2, +1);
         if (k3 <= il_MaxK) RemoveTermsSmallPrime(p3, k3, +1);
         if (k4 <= il_MaxK) RemoveTermsSmallPrime(p4, k4, +1);
      }
      else
      {
         // This function allows for a single term to be removed
         if (k1 <= il_MaxK) RemoveTermsBigPrime(p1, k1, +1);
         if (k2 <= il_MaxK) RemoveTermsBigPrime(p2, k2, +1);
         if (k3 <= il_MaxK) RemoveTermsBigPrime(p3, k3, +1);
         if (k4 <= il_MaxK) RemoveTermsBigPrime(p4, k4, +1);
      }

      SetLargestPrimeTested(p4, 4);
   }
         
   // Adjust for the possibility that we tested the same prime
   // more than once at the end of the list.
   if (p4 == p3)
      SetLargestPrimeTested(p4, -1);
   if (p4 == p2)
      SetLargestPrimeTested(p4, -1);
   if (p4 == p1)
      SetLargestPrimeTested(p4, -1);
}

void  TwinWorker::TestLargeB(void)
{
   uint64_t ks[4], ps[4];
   uint32_t idx;
         
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
   
   
      if (it_TermType == TT_BN)
      {
         ks[0] = InvMod32(ii_Base, ps[0]);
         ks[1] = InvMod32(ii_Base, ps[1]);
         ks[2] = InvMod32(ii_Base, ps[2]);
         ks[3] = InvMod32(ii_Base, ps[3]);
         
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

      if (ps[0] <= il_MaxK)
      {
         // This function allows for multiple terms to be removed
         if (ks[0] <= il_MaxK) RemoveTermsSmallPrime(ps[0], ks[0], -1);
         if (ks[1] <= il_MaxK) RemoveTermsSmallPrime(ps[1], ks[1], -1);
         if (ks[2] <= il_MaxK) RemoveTermsSmallPrime(ps[2], ks[2], -1);
         if (ks[3] <= il_MaxK) RemoveTermsSmallPrime(ps[3], ks[3], -1);
      }
      else
      {
         // This function allows for a single term to be removed
         if (ks[0] <= il_MaxK) RemoveTermsBigPrime(ps[0], ks[0], -1);
         if (ks[1] <= il_MaxK) RemoveTermsBigPrime(ps[1], ks[1], -1);
         if (ks[2] <= il_MaxK) RemoveTermsBigPrime(ps[2], ks[2], -1);
         if (ks[3] <= il_MaxK) RemoveTermsBigPrime(ps[3], ks[3], -1);
      }

      ks[0] = ps[0] - ks[0];
      ks[1] = ps[1] - ks[1];
      ks[2] = ps[2] - ks[2];
      ks[3] = ps[3] - ks[3];
      
      if (ps[0] <= il_MaxK)
      {
         // This function allows for multiple terms to be removed
         if (ks[0] <= il_MaxK) RemoveTermsSmallPrime(ps[0], ks[0], +1);
         if (ks[1] <= il_MaxK) RemoveTermsSmallPrime(ps[1], ks[1], +1);
         if (ks[2] <= il_MaxK) RemoveTermsSmallPrime(ps[2], ks[2], +1);
         if (ks[3] <= il_MaxK) RemoveTermsSmallPrime(ps[3], ks[3], +1);
      }
      else
      {
         // This function allows for a single term to be removed
         if (ks[0] <= il_MaxK) RemoveTermsBigPrime(ps[0], ks[0], +1);
         if (ks[1] <= il_MaxK) RemoveTermsBigPrime(ps[1], ks[1], +1);
         if (ks[2] <= il_MaxK) RemoveTermsBigPrime(ps[2], ks[2], +1);
         if (ks[3] <= il_MaxK) RemoveTermsBigPrime(ps[3], ks[3], +1);
      }

      SetLargestPrimeTested(ps[3], 4);
   }
}

// This must be used when prime <= il_MaxK.
void    TwinWorker::RemoveTermsSmallPrime(uint64_t prime, uint64_t k, int32_t c)
{
   // This primes will yield no factor because b^n = 0 (mod p)
   if (ii_Base % prime == 0)
      return;
   
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
   
   if (k > il_MaxK)
      return;
   
   ip_TwinApp->ReportFactor(prime, k, c);
}


// Using this bypasses a number of if checks that can be done when prime > il_MaxK.
void    TwinWorker::RemoveTermsBigPrime(uint64_t prime, uint64_t k, int32_t c)
{   
   // Make sure that k >= il_MinK
   if (k < il_MinK)
      return;
   
   ip_TwinApp->ReportFactor(prime, k, c);
}

void    TwinWorker::BuildBaseInverses(void)
{
	ii_BaseInverses = (uint32_t *) xmalloc(ii_Base * sizeof(uint32_t));

	ii_BaseInverses[0] = 0;
	ii_BaseInverses[1] = ii_Base-1;

	for (uint32_t i = 2; i < ii_Base; i++)
		ii_BaseInverses[i] = EuclidExtendedGCD(i, ii_Base);
}

uint32_t TwinWorker::EuclidExtendedGCD(uint32_t a, uint32_t base)
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

// Determine the minimum prime between kmin*b^n+c and kmax*b^n+c so
// that we will avoid removing those terms in that range that are prime
void  TwinWorker::DeterminePrimeTermRange(void)
{  
   uint64_t  maxPrime = ip_TwinApp->GetMaxPrime();
   double    tooBig = (double) KMAX_MAX;
   double    b = (double) ii_Base;
   double    bpown = 1.0;
   double    mink = (double) il_MinK;
   
   // KMAX_MAX is the same as PMAX_MAX, 2^62.
   if (maxPrime == 0)
      tooBig = (double) KMAX_MAX;
   else
      tooBig = (double) maxPrime;
   
   il_BpowN = 0;
   
   for (uint32_t i=0; i<ii_N; i++)
   {
      bpown *= b;

      // If mink*b^n-1.0 > 2^62, then we don't need to worry about removing terms
      // that are prime.
      if (mink*bpown-1.0 > tooBig)
         return;
   }

   il_BpowN = (uint64_t) bpown;
}
