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
   vector<uint64_t>::iterator it = iv_Primes.begin();
   uint64_t p0 = *it;
   
   if (p0 < il_MaxK || p0 < 50)
      TestMegaPrimeChunkSmall();
   else
      TestMegaPrimeChunkLarge();
}

void  GFNDivisorWorker::TestMegaPrimeChunkSmall(void)
{
   uint64_t k1, k2, k3, k4, ks[4];
   uint64_t p1, p2, p3, p4, ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t n1, n2, n3, n4;
   uint32_t bits1, bits2, bits3, bits4;
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
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
      
      n1 = n2 = n3 = n4 = ii_MinN;

      while (n1 <= ii_MaxN)
      {
         bits1 = __builtin_ctzll(k1);
         bits2 = __builtin_ctzll(k2);
         bits3 = __builtin_ctzll(k3);
         bits4 = __builtin_ctzll(k4);
         
         k1 >>= bits1;
         k2 >>= bits2;
         k3 >>= bits3;
         k4 >>= bits4;
         
         n1 += bits1;
         n2 += bits2;
         n3 += bits3;
         n4 += bits4;
         
         if (k1 >= il_MinK && k1 <= il_MaxK && n1 <= ii_MaxN) RemoveTermsSmallPrime(k1, n1, p1);
         if (k2 >= il_MinK && k2 <= il_MaxK && n2 <= ii_MaxN) RemoveTermsSmallPrime(k2, n2, p2);
         if (k3 >= il_MinK && k3 <= il_MaxK && n3 <= ii_MaxN) RemoveTermsSmallPrime(k3, n3, p3);
         if (k4 >= il_MinK && k4 <= il_MaxK && n4 <= ii_MaxN) RemoveTermsSmallPrime(k4, n4, p4);
         
         k1 += p1;
         k2 += p2;
         k3 += p3;
         k4 += p4;
      }

      while (n2 <= ii_MaxN)
      {
         bits2 = __builtin_ctzll(k2);
         
         k2 >>= bits2;
         n2 += bits2;
         
         if (k2 >= il_MinK && k2 <= il_MaxK) RemoveTermsSmallPrime(k2, n2, p2);
         
         k2 += p2;
      }

      while (n3 <= ii_MaxN)
      {
         bits3 = __builtin_ctzll(k3);
         
         k3 >>= bits3;
         n3 += bits3;
         
         if (k3 >= il_MinK && k3 <= il_MaxK && n3 <= ii_MaxN) RemoveTermsSmallPrime(k3, n3, p3);
         
         k3 += p3;
      }   

      while (n4 <= ii_MaxN)
      {
         bits4 = __builtin_ctzll(k4);

         k4 >>= bits4;
         n4 += bits4;
                  
         if (k4 >= il_MinK && k4 <= il_MaxK && n4 <= ii_MaxN) RemoveTermsSmallPrime(k4, n4, p4);
         
         k4 += p4;
      }   

      SetLargestPrimeTested(p4, 4);
   
      if (p4 >= maxPrime)
         break;
   }
}

void  GFNDivisorWorker::TestMegaPrimeChunkLarge(void)
{
   uint64_t k1, k2, k3, k4, ks[4];
   uint64_t p1, p2, p3, p4, ps[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t n1, n2, n3, n4;
   uint32_t bits1, bits2, bits3, bits4;
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
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
      
      n1 = n2 = n3 = n4 = ii_MinN;

      while (n1 <= ii_MaxN)
      {
         bits1 = __builtin_ctzll(k1);
         bits2 = __builtin_ctzll(k2);
         bits3 = __builtin_ctzll(k3);
         bits4 = __builtin_ctzll(k4);
         
         k1 >>= bits1;
         k2 >>= bits2;
         k3 >>= bits3;
         k4 >>= bits4;
         
         n1 += bits1;
         n2 += bits2;
         n3 += bits3;
         n4 += bits4;
          
         if (k1 >= il_MinK && k1 <= il_MaxK && n1 <= ii_MaxN) RemoveTermsBigPrime(k1, n1, p1);
         if (k2 >= il_MinK && k2 <= il_MaxK && n2 <= ii_MaxN) RemoveTermsBigPrime(k2, n2, p2);
         if (k3 >= il_MinK && k3 <= il_MaxK && n3 <= ii_MaxN) RemoveTermsBigPrime(k3, n3, p3);
         if (k4 >= il_MinK && k4 <= il_MaxK && n4 <= ii_MaxN) RemoveTermsBigPrime(k4, n4, p4);
         
         k1 += p1;
         k2 += p2;
         k3 += p3;
         k4 += p4;
      }

      while (n2 <= ii_MaxN)
      {
         bits2 = __builtin_ctzll(k2);
         
         k2 >>= bits2;
         n2 += bits2;
         
         if (k2 >= il_MinK && k2 <= il_MaxK) RemoveTermsBigPrime(k2, n2, p2);
         
         k2 += p2;
      }

      while (n3 <= ii_MaxN)
      {
         bits3 = __builtin_ctzll(k3);
         
         k3 >>= bits3;
         n3 += bits3;
         
         if (k3 >= il_MinK && k3 <= il_MaxK) RemoveTermsBigPrime(k3, n3, p3);
         
         k3 += p3;
      }   

      while (n4 <= ii_MaxN)
      {
         bits4 = __builtin_ctzll(k4);

         k4 >>= bits4;
         n4 += bits4;
                  
         if (k4 >= il_MinK && k4 <= il_MaxK) RemoveTermsBigPrime(k4, n4, p4);
         
         k4 += p4;
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
// Do not report k/n combinations if the k*2^n+1 is divisible by any p < 50
void    GFNDivisorWorker::RemoveTermsBigPrime(uint64_t k, uint32_t n, uint64_t prime)
{   
   uint32_t smallN;
   uint64_t smallK;

   smallN = n % (2);
   smallK = k % (3);
   if ((smallK << smallN) % (3) == 2) return;
     
   smallN = n % (4);
   smallK = k % (5);
   if ((smallK << smallN) % (5) == 4) return;
   
   smallN = n % (6);
   smallK = k % (7);
   if ((smallK << smallN) % (7) == 6) return;
   
   smallN = n % (10);
   smallK = k % (11);
   if ((smallK << smallN) % (11) == 10) return;
   
   smallN = n % (12);
   smallK = k % (13);
   if ((smallK << smallN) % (13) == 12) return;
   
   smallN = n % (16);
   smallK = k % (17);
   if ((smallK << smallN) % (17) == 16) return;
   
   smallN = n % (18);
   smallK = k % (19);
   if ((smallK << smallN) % (19) == 18) return;
   
   smallN = n % (22);
   smallK = k % (23);
   if ((smallK << smallN) % (23) == 22) return;
   
   smallN = n % (28);
   smallK = k % (29);
   if ((smallK << smallN) % (29) == 28) return;
   
   smallN = n % (30);
   smallK = k % (31);
   if ((smallK << smallN) % (31) == 30) return;

   smallN = n % (36);
   smallK = k % (37);
   if ((smallK << smallN) % (37) == 36) return;
   
   smallN = n % (40);
   smallK = k % (41);
   if ((smallK << smallN) % (41) == 40) return;
   
   smallN = n % (42);
   smallK = k % (43);
   if ((smallK << smallN) % (43) == 42) return;
   
   smallN = n % (46);
   smallK = k % (47);
   if ((smallK << smallN) % (47) == 46) return;

   ip_GFNDivisorApp->ReportFactor(prime, k, n, true);
}
