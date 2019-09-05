/* DMDivisorWorker.cpp -- (C) Mark Rodenkirch, September 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "DMDivisorWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

DMDivisorWorker::DMDivisorWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_DMDivisorApp = (DMDivisorApp *) theApp;
   
   il_MinK = ip_DMDivisorApp->GetMinK();
   il_MaxK = ip_DMDivisorApp->GetMaxK();
   ii_Exp = ip_DMDivisorApp->GetN();
   
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
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   while (it != iv_Primes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;
      
      bs[0] = bs[1] = bs[2] = bs[3] = 2;
      
      fpu_powmod_4b_1n_4p(bs, ii_Exp, ps);

      // Now bs = 2^exp-1 (mod p)
      bs[0]--;
      bs[1]--;
      bs[2]--;
      bs[3]--;
      
      // We are looking for p such that 2*k*bs+1 (mod p) = 0
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
      
      // We have now solved for k
         
      if (ps[0] <= il_MaxK)
      {
         if (k1 <= il_MaxK) RemoveTermsSmallPrime(ps[0], k1);
         if (k2 <= il_MaxK) RemoveTermsSmallPrime(ps[1], k2);
         if (k3 <= il_MaxK) RemoveTermsSmallPrime(ps[2], k3);
         if (k4 <= il_MaxK) RemoveTermsSmallPrime(ps[3], k4);
      }
      else
      {
         if (k1 <= il_MaxK) RemoveTermsBigPrime(ps[0], k1);
         if (k2 <= il_MaxK) RemoveTermsBigPrime(ps[1], k2);
         if (k3 <= il_MaxK) RemoveTermsBigPrime(ps[2], k3);
         if (k4 <= il_MaxK) RemoveTermsBigPrime(ps[3], k4);
      }

      SetLargestPrimeTested(ps[3], 4);
   
      if (ps[3] >= maxPrime)
         break;
   }
}

void  DMDivisorWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("DMDivisorWorker::TestMiniPrimeChunk not implemented");
}

// This must be used when prime <= il_MaxK.
void    DMDivisorWorker::RemoveTermsSmallPrime(uint64_t prime, uint64_t k)
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
   
   if (k > il_MaxK)
      return;
   
   do
	{
      if (ip_DMDivisorApp->ReportFactor(prime, k))
         VerifyFactor(prime, k);

		k += prime; 
	} while (k <= il_MaxK);
   
}

// Using this bypasses a number of if checks that can be done when prime > il_MaxK.
void    DMDivisorWorker::RemoveTermsBigPrime(uint64_t prime, uint64_t k)
{         
   // Make sure that k >= il_MinK
   if (k < il_MinK)
      return;
   
   if (ip_DMDivisorApp->ReportFactor(prime, k))
      VerifyFactor(prime, k);
}

void  DMDivisorWorker::VerifyFactor(uint64_t prime, uint64_t k)
{
   uint64_t rem;

   fpu_push_1divp(prime);
      
   rem = fpu_powmod(2, ii_Exp, prime);
   
   if (rem == 0)
      rem = prime - 1;
   else
      rem--;
   
   rem = fpu_mulmod(rem, k, prime);

   rem <<= 1;
   rem++;
   
   if (rem >= prime)
      rem -= prime;
   
   if (rem != 0)
      FatalError("2*%" PRIu64"*(2^%" PRIu64"-1)+1 mod %" PRIu64" = %" PRIu64"", k, ii_Exp, prime, rem);
   
   fpu_pop();
}
