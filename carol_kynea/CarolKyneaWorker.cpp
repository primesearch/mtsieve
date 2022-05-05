/* CarolKynea.cpp -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include "CarolKyneaWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

CarolKyneaWorker::CarolKyneaWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_CarolKyneaApp = (CarolKyneaApp *) theApp;
   
   ii_Base = ip_CarolKyneaApp->GetBase();
   ii_MinN = ip_CarolKyneaApp->GetMinN();
   ii_MaxN = ip_CarolKyneaApp->GetMaxN();

   uint32_t r = ii_MaxN - ii_MinN + 1;

   // In the worst case we will do do one table insertion and one mulmod
   // for m baby steps, then s table lookups and s mulmods for M giant
   // steps. The average case depends on how many solutions are found
   // and how early in the loop they are found, which I don't know how
   // to analyse. However for the worst case we just want to minimise
   // m + s*M subject to m*M >= r, which is when m = sqrt(s*r).
  
   ii_GiantSteps = MAX(1, sqrt((double) r/ROOT_COUNT));
   ii_BabySteps = MIN(r, ceil((double) r/ii_GiantSteps));

   if (ii_BabySteps > HASH_MAX_ELTS)
   {
      ii_GiantSteps = ceil((double)r/HASH_MAX_ELTS);
      ii_BabySteps = ceil((double)r/ii_GiantSteps);
   }

   ii_SieveLow = ip_CarolKyneaApp->GetMinN();
   ii_SieveRange = ii_BabySteps*ii_GiantSteps;
   
   assert(ii_SieveLow <= ip_CarolKyneaApp->GetMinN());
   assert(ip_CarolKyneaApp->GetMaxN() < ii_SieveLow+ii_SieveRange);
   
   ip_HashTable = new HashTable(ii_BabySteps);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CarolKyneaWorker::CleanUp(void)
{
   delete ip_HashTable;
}

void  CarolKyneaWorker::TestMegaPrimeChunk(void)
{
   uint64_t  maxPrime = ip_App->GetMaxPrime();
   uint64_t  thePrime = 0, root1, root2;
   vector<uint64_t>::iterator it = iv_Primes.begin();

   while (it != iv_Primes.end())
   {
      thePrime = *it;
      it++;
      
      if (ii_Base % thePrime == 0)
         continue;
      
      // Skip this prime if there are no values x such that x^2 = 2 (mod p)
      if (!IsQuadraticResidue(2, thePrime))
         continue;
   
      fpu_push_1divp(thePrime);

      // Find root of x^2 = 2 (mod p)
      root1 = FindRoot(thePrime);

      root2 = thePrime - root1;

      // It is possible that findRoot returns where x^2 = -2 (mod p)
      if (fpu_mulmod(root1, root1, thePrime) != 2)         
         ip_CarolKyneaApp->WriteToConsole(COT_SIEVE, "%" PRIu64" is not a root (mod %" PRIu64")", root1, thePrime);
      
      if (fpu_mulmod(root2, root2, thePrime) != 2) 
         ip_CarolKyneaApp->WriteToConsole(COT_SIEVE, "%" PRIu64" is not a root (mod %" PRIu64")", root2, thePrime);
    
      io_Sequence[0].root = root1 - 1;
      io_Sequence[0].c    = +1;
      io_Sequence[1].root = root2 - 1;
      io_Sequence[1].c    = +1;
      io_Sequence[2].root = root1 + 1;
      io_Sequence[2].c    = -1;
      io_Sequence[3].root = root2 + 1;
      io_Sequence[3].c    = -1;
         
      DiscreteLog(thePrime);
      
      fpu_pop();

      SetLargestPrimeTested(thePrime, 1);
      
      if (thePrime >= maxPrime)
         break;
   }
}

void  CarolKyneaWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CarolKyneaWorker::TestMiniPrimeChunk not implemented");
}

void  CarolKyneaWorker::DiscreteLog(uint64_t p)
{
   uint64_t  b, bj0;
   uint32_t  i, j, k;
   uint64_t  inv_pb;

   b = ii_Base % p;
   
   ip_HashTable->Clear();
   
   // Precompute 1/b^d (mod p) for 0 <= d <= Q.
   inv_pb = InvMod32(b, p);
   
   if (inv_pb == 0)
      return;

   for (i = 0; i < ROOT_COUNT; i++)
      il_A[i] = io_Sequence[i].root % p;

   b = ii_Base;
   bj0 = fpu_powmod(b, ii_MinN, p);
   
   i = BabySteps(b, bj0, p);
   
   if (i > 0)
   {
      // i is the order of b (mod p). This is all the information we need to
      // determine every solution for this p, no giant steps are needed.
      for (k = 0; k < ROOT_COUNT; k++)
         for (j = ip_HashTable->Lookup(il_A[k]); j < ii_SieveRange; j += i)
            CheckFactor(p, ii_SieveLow+j, io_Sequence[k].c);
         
      return;
   }
   
   // First giant step
   for (k = 0; k < ROOT_COUNT; k++)
      if ((j = ip_HashTable->Lookup(il_A[k])) != HASH_NOT_FOUND)
        CheckFactor(p, ii_SieveLow+j, io_Sequence[k].c);

   // Remaining giant steps
   b = fpu_powmod(inv_pb, ii_BabySteps, p); /* b <- 1/b^m (mod p) */

   fpu_push_adivb(b, p);
         
   for (i = 1; i < ii_GiantSteps; i++)
   {
      fpu_mulmod_iter_4a(il_A, b, p);

      if ((j = ip_HashTable->Lookup(il_A[0])) != HASH_NOT_FOUND)
         CheckFactor(p, ii_SieveLow+i*ii_BabySteps+j, io_Sequence[0].c);
         
      if ((j = ip_HashTable->Lookup(il_A[1])) != HASH_NOT_FOUND)
         CheckFactor(p, ii_SieveLow+i*ii_BabySteps+j, io_Sequence[1].c);
      
      if ((j = ip_HashTable->Lookup(il_A[2])) != HASH_NOT_FOUND)
         CheckFactor(p, ii_SieveLow+i*ii_BabySteps+j, io_Sequence[2].c);
      
      if ((j = ip_HashTable->Lookup(il_A[3])) != HASH_NOT_FOUND)
         CheckFactor(p, ii_SieveLow+i*ii_BabySteps+j, io_Sequence[3].c);
   }

   fpu_pop();
}

uint32_t  CarolKyneaWorker::BabySteps(uint64_t b, uint64_t bj0, uint64_t p)
{
   uint64_t bj;
   uint32_t j;
   
   fpu_push_adivb(b, p);
  
   for (j = 0, bj = bj0; j < ii_BabySteps; j++)
   {
      ip_HashTable->Insert(bj, j);
      
      bj = fpu_mulmod_iter(bj, b, p);
            
      if (bj == bj0)
      {
         fpu_pop();
         return j+1;
      }
   }

   fpu_pop();
   return 0;
}

// Find x such that x^2 = 2 mod p
// This is solved using Hensel's Lemma
uint64_t	CarolKyneaWorker::FindRoot(uint64_t p)
{
	uint64_t	   i, s, t, d, m, rem, A, D;

	if ((p & 7) == 7)
		return fpu_powmod(2, (p+1) >> 2, p);

	t = p - 1;
	s = 0;
	while (!(t & 1))
	{
		s++;
		t >>= 1;
	}

	A = fpu_powmod(2, t, p);

	// Find value d where Lengendre Symbol is -1
	for (d=3; d<p; d++)
		if (!IsQuadraticResidue(d, p))
			break;

	D = fpu_powmod(d, t, p);

	m = 0;
	for (i=0; i<s; i++)
	{
		if (m == 0)
			rem = 1;
		else
			rem = fpu_powmod(D, m, p);

		rem = fpu_mulmod(rem, A, p);
		rem = fpu_powmod(rem, 1 << (s - 1 - i), p);
		
		if (rem == p - 1)
			m += (1 << i);
	}

	if (m == 0)
		rem = 1;
	else
		rem = fpu_powmod(D, m >> 1, p);

	i = fpu_powmod(2, (t+1) >> 1, p);
   
	return fpu_mulmod(rem, i, p);
}

void CarolKyneaWorker::CheckFactor(uint64_t p, uint32_t n, int32_t c)
{
   if (n < 62)
   {
      uint64_t term = 1;
      uint64_t maxBeforeOverflow = (1 << 31);
      uint32_t i = 0;
      bool isFactor = false;
      
      do
      {
         term *= ii_Base;
         i++;
         
         if (term + c > maxBeforeOverflow)
            isFactor = true;
         
         if (((term + c) * (term + c)) > p + 2)
            isFactor = true;
      } while (i < n && !isFactor);

      if (!isFactor)
      {
         ip_CarolKyneaApp->ReportPrime(p, n, c);
         return;
      }
   }
   
   if (ip_CarolKyneaApp->ReportFactor(p, n, c))
      VerifyFactor(p, n, c);
}

void CarolKyneaWorker::VerifyFactor(uint64_t p, uint32_t n, int32_t c)
{
   uint64_t rem;
   
   fpu_push_1divp(p);
   
   rem = fpu_powmod(ii_Base, n, p);
   
   rem = fpu_mulmod(rem + c, rem + c, p);
   
   fpu_pop();
   
   if (rem != 2)
      FatalError("(%u^%u%+d)-2 mod %" PRIu64" = %" PRIu64"", ii_Base, n, c, p, rem-2);
}