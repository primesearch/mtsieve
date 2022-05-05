/* SophieGermainWorker.cpp -- (C) Mark Rodenkirch, July 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "SophieGermainWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

SophieGermainWorker::SophieGermainWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_SophieGermainApp = (SophieGermainApp *) theApp;
   
   il_MinK = ip_SophieGermainApp->GetMinK();
   il_MaxK = ip_SophieGermainApp->GetMaxK();
   ii_N = ip_SophieGermainApp->GetN();
     
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  SophieGermainWorker::CleanUp(void)
{
}

void  SophieGermainWorker::TestMegaPrimeChunk(void)
{
   uint64_t k1, k2, k3, k4, ks[4];
   uint64_t p1, p2, p3, p4 = 0, ps[4];
   uint64_t p0, maxPrime = ip_App->GetMaxPrime();
   
   vector<uint64_t>::iterator it = iv_Primes.begin();

   p0 = *it;
         
   while (it != iv_Primes.end())
   {
      ps[0] = p1 = *it;
      ks[0] = k1 = (1+p1) >> 1;
      it++;
      
      ps[1] = p2 = *it;
      ks[1] = k2 = (1+p2) >> 1;
      it++;
      
      ps[2] = p3 = *it;
      ks[2] = k3 = (1+p3) >> 1;
      it++;
      
      ps[3] = p4 = *it;
      ks[3] = k4 = (1+p4) >> 1;
      it++;      
      
      // Starting with k*2^n = 1 (mod p) 
      //           --> k = (1/2)^n (mod p)
      //           --> k = inverse^n (mod p)
      fpu_powmod_4b_1n_4p(ks, ii_N, ps);
      
      k1 = p1 - ks[0];
      k2 = p2 - ks[1];
      k3 = p3 - ks[2];
      k4 = p4 - ks[3];

      if (p0 <= il_MaxK)
      {
         if (k1 <= il_MaxK) RemoveTermsSmallPrime(p1, k1, true);
         if (k2 <= il_MaxK) RemoveTermsSmallPrime(p2, k2, true);
         if (k3 <= il_MaxK) RemoveTermsSmallPrime(p3, k3, true);
         if (k4 <= il_MaxK) RemoveTermsSmallPrime(p4, k4, true);
      }
      else
      {
         if (k1 <= il_MaxK) RemoveTermsLargePrime(p1, k1, true);
         if (k2 <= il_MaxK) RemoveTermsLargePrime(p2, k2, true);
         if (k3 <= il_MaxK) RemoveTermsLargePrime(p3, k3, true);
         if (k4 <= il_MaxK) RemoveTermsLargePrime(p4, k4, true);
      }

      // For the next n, divide k by 2
      // Note that k*2^n+1 = (k/2)*2^(n+1)+1
      
      // Make sure k is even first
      if (k1 & 1) k1 += p1;
      if (k2 & 1) k2 += p2;
      if (k3 & 1) k3 += p3;
      if (k4 & 1) k4 += p4;
         
      k1 >>= 1;
      k2 >>= 1;
      k3 >>= 1;
      k4 >>= 1;
      
      if (p0 <= il_MaxK)
      {
         if (k1 <= il_MaxK) RemoveTermsSmallPrime(p1, k1, false);
         if (k2 <= il_MaxK) RemoveTermsSmallPrime(p2, k2, false);
         if (k3 <= il_MaxK) RemoveTermsSmallPrime(p3, k3, false);
         if (k4 <= il_MaxK) RemoveTermsSmallPrime(p4, k4, false);
      }
      else
      {
         if (k1 <= il_MaxK) RemoveTermsLargePrime(p1, k1, false);
         if (k2 <= il_MaxK) RemoveTermsLargePrime(p2, k2, false);
         if (k3 <= il_MaxK) RemoveTermsLargePrime(p3, k3, false);
         if (k4 <= il_MaxK) RemoveTermsLargePrime(p4, k4, false);
      }


      SetLargestPrimeTested(p4, 4);
   
      if (p4 >= maxPrime)
         break;
   }
}

void  SophieGermainWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SophieGermainWorker::TestMiniPrimeChunk not implemented");
}

// This must be used when prime <= il_MaxK.
void    SophieGermainWorker::RemoveTermsSmallPrime(uint64_t prime, uint64_t k, bool firstOfPair)
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

   // Make sure that k is odd
   if (!(k & 1))
      k += prime;

   if (k > il_MaxK)
      return;
      
   fpu_push_1divp(prime);
   il_BpowN = fpu_powmod(2, ii_N, prime);
   
   do
	{
      if (ip_SophieGermainApp->ReportFactor(prime, k, firstOfPair))
         VerifyFactor(prime, k, firstOfPair);

		k += (prime << 1); 
	} while (k <= il_MaxK);
   
   fpu_pop();
}

// Using this bypasses a number of if checks that can be done when prime > il_MaxK.
void    SophieGermainWorker::RemoveTermsLargePrime(uint64_t prime, uint64_t k, bool firstOfPair)
{
   // Make sure that k >= il_MinK
   if (k < il_MinK)
      return;

   // Make sure that k is odd
   if (!(k & 1))
      return;

   if (ip_SophieGermainApp->ReportFactor(prime, k, firstOfPair))
      VerifyExternalFactor(prime, k, ii_N, firstOfPair);
}

void  SophieGermainWorker::VerifyExternalFactor(uint64_t prime, uint64_t k, uint32_t n, bool firstOfPair)
{
   fpu_push_1divp(prime);
   
   il_BpowN = fpu_powmod(2, ii_N, prime);
   
   VerifyFactor(prime, k, firstOfPair);
   
   fpu_pop();
}

void  SophieGermainWorker::VerifyFactor(uint64_t prime, uint64_t k, bool firstOfPair)
{
   uint64_t rem;

   rem = fpu_mulmod(il_BpowN, k, prime);

   if (!firstOfPair)
   {
      rem <<= 1;

      if (rem >= prime)
         rem -= prime;
   }
   
   if (rem != prime - 1)
   {
      if (firstOfPair)
         FatalError("%" PRIu64"*2^%u+1 mod %" PRIu64" = %" PRIu64"", k, ii_N, prime, rem + 1);
      else
         FatalError("2*(%" PRIu64"*2^%u+1)-1 mod %" PRIu64" = %" PRIu64"", k, ii_N+1, prime, rem + 1);
   }
}