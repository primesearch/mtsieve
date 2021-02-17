/* GFNDivisorWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>

#include "GFNDivisorWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

GFNDivisorWorker::GFNDivisorWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_GFNDivisorApp = (GFNDivisorApp *) theApp;
   
   il_KCount = ip_GFNDivisorApp->GetKCount();
   ii_NCount = ip_GFNDivisorApp->GetNCount();
   
   il_MinK = ip_GFNDivisorApp->GetMinK();
   il_MaxK = ip_GFNDivisorApp->GetMaxK();
   
   ii_MinN = ip_GFNDivisorApp->GetMinN();
   ii_MaxN = ip_GFNDivisorApp->GetMaxN();
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  GFNDivisorWorker::CleanUp(void)
{
}

void  GFNDivisorWorker::TestMegaPrimeChunk(void)
{
   uint64_t k1, k2, k3, k4, ks[4];
   uint64_t p1, p2, p3, p4 = 0, ps[4];
   uint64_t p0, maxPrime = ip_App->GetMaxPrime();
   uint32_t n;
   
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
      fpu_powmod_4b_1n_4p(ks, ii_MinN, ps);
      
      k1 = p1 - ks[0];
      k2 = p2 - ks[1];
      k3 = p3 - ks[2];
      k4 = p4 - ks[3];

      // look at app_thread_fun.c for speed ups.
      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (p0 <= il_MaxK)
         {
            if (k1 <= il_MaxK) RemoveTermsSmallPrime(k1, n, p1);
            if (k2 <= il_MaxK) RemoveTermsSmallPrime(k2, n, p2);
            if (k3 <= il_MaxK) RemoveTermsSmallPrime(k3, n, p3);
            if (k4 <= il_MaxK) RemoveTermsSmallPrime(k4, n, p4);
         }
         else
         {
            if (k1 <= il_MaxK) RemoveTermsBigPrime(k1, n, p1);
            if (k2 <= il_MaxK) RemoveTermsBigPrime(k2, n, p2);
            if (k3 <= il_MaxK) RemoveTermsBigPrime(k3, n, p3);
            if (k4 <= il_MaxK) RemoveTermsBigPrime(k4, n, p4);
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
      }

      SetLargestPrimeTested(p4, 4);
   
      if (p4 >= maxPrime)
         break;
   }
}

void  GFNDivisorWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GFNDivisorWorker::TestMiniPrimeChunk not implemented");
}

// This must be used when prime <= il_MaxK.
void    GFNDivisorWorker::RemoveTermsSmallPrime(uint64_t k, uint32_t n, uint64_t prime)
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

   uint32_t verifiedCount = 0;
   
   do
	{
      verifiedCount++;
      
      // If the first few are valid, then assume that the rest are valid.  This will
      // speed up testing of small primes.
      ip_GFNDivisorApp->ReportFactor(prime, k, n, (verifiedCount < 5));

		k += (prime << 1); 
	} while (k <= il_MaxK);
}

// Using this bypasses a number of if checks that can be done when prime > il_MaxK.
void    GFNDivisorWorker::RemoveTermsBigPrime(uint64_t k, uint32_t n, uint64_t prime)
{         
   // Make sure that k >= il_MinK
   if (k < il_MinK)
      return;

   // Make sure that k is odd
   if (!(k & 1))
      return;
   
   ip_GFNDivisorApp->ReportFactor(prime, k, n, true);
}
