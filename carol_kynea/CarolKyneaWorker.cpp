/* CarolKynea.cpp -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include "CarolKyneaWorker.h"

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
   uint64_t  thePrime = 0;
   uint64_t  inv_pb;

   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx++)
   {
      thePrime = il_PrimeList[pIdx+0];
      
      if (ii_Base % thePrime == 0)
         continue;
      
      // Skip this prime if there are no values x such that x^2 = 2 (mod p)
      if (!IsQuadraticResidue(2, thePrime))
         continue;
   
      // Precompute 1/b^d (mod p) for 0 <= d <= Q.
      inv_pb = InvMod32(ii_Base % thePrime, thePrime);
      
      if (inv_pb == 0)
         return;
      
      MpArith mp(thePrime);
      MpRes   res2 = mp.nToRes(2);
   
      // Find root of x^2 = 2 (mod p)
      ir_Root1 = FindRoot(thePrime, mp);

      ir_Root2 = mp.sub(mp.nToRes(thePrime), ir_Root1); 

      // It is possible that findRoot returns where x^2 = -2 (mod p)
      if (mp.mul(ir_Root1, ir_Root1) != res2)
         ip_CarolKyneaApp->WriteToConsole(COT_SIEVE, "%" PRIu64" is not a root (mod %" PRIu64")", mp.resToN(ir_Root1), thePrime);
      
      if (mp.mul(ir_Root2, ir_Root2) != res2)
         ip_CarolKyneaApp->WriteToConsole(COT_SIEVE, "%" PRIu64" is not a root (mod %" PRIu64")", mp.resToN(ir_Root2), thePrime);

      DiscreteLog(thePrime, mp, inv_pb);
      
      SetLargestPrimeTested(thePrime, 1);
      
      if (thePrime >= maxPrime)
         break;
   }
}

void  CarolKyneaWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CarolKyneaWorker::TestMiniPrimeChunk not implemented");
}

// Find x such that x^2 = 2 mod p
// This is solved using Hensel's Lemma
uint64_t	CarolKyneaWorker::FindRoot(uint64_t p, MpArith mp)
{
	uint64_t	   i, s, t, d, m;
   MpRes       res, resA, resD;
   MpRes       res2 = mp.nToRes(2);
   MpRes       resPM1 = mp.nToRes(p-1);
   
	if ((p & 7) == 7)
      return mp.pow(res2, (p+1) >> 2);
   
	t = p - 1;
	s = 0;
	while (!(t & 1))
	{
		s++;
		t >>= 1;
	}

   resA = mp.pow(res2, t);

	// Find value d where Lengendre Symbol is -1
	for (d=3; d<p; d++)
		if (!IsQuadraticResidue(d, p))
			break;

   resD = mp.nToRes(d);
   
   resD = mp.pow(resD, t);

	m = 0;
	for (i=0; i<s; i++)
	{
		if (m == 0)
			res = mp.one();
		else
			res = mp.pow(resD, m);

      res = mp.mul(res, resA);
      res = mp.pow(res, 1 << (s - 1 - i));
		
		if (res == resPM1)
			m += (1 << i);
	}

	if (m == 0)
		res = mp.one();
	else
		res = mp.pow(resD, m >> 1);
   
   resA = mp.pow(res2, (t+1) >> 1);
   
   return mp.mul(res, resA);
}

void  CarolKyneaWorker::DiscreteLog(uint64_t p, MpArith mp, uint64_t inv_pb)
{
   uint32_t  orderOfB;
   MpRes     resBase = mp.nToRes(ii_Base);
   MpRes     resA[ROOT_COUNT];
   uint32_t  jHash[ROOT_COUNT];
   
   ip_HashTable->Clear();
      
   if (inv_pb == 0)
      return;

   resA[0] = mp.sub(ir_Root1, mp.one());
   resA[1] = mp.sub(ir_Root2, mp.one());
   resA[2] = mp.add(ir_Root1, mp.one());
   resA[3] = mp.add(ir_Root2, mp.one());

   orderOfB = BabySteps(p, mp, resBase);

   jHash[0] = ip_HashTable->Lookup(resA[0]);
   jHash[1] = ip_HashTable->Lookup(resA[1]);
   jHash[2] = ip_HashTable->Lookup(resA[2]);
   jHash[3] = ip_HashTable->Lookup(resA[3]);
   
   if (orderOfB > 0)
   {
      uint32_t j;
      
      // orderOfB is the order of b (mod p). This is all the information we need to
      // determine every solution for this p, no giant steps are needed.
      
      for (j = jHash[0]; j < ii_SieveRange; j += orderOfB)
         ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+j, +1);
      
      for (j = jHash[1]; j < ii_SieveRange; j += orderOfB)
         ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+j, +1);

      for (j = jHash[2]; j < ii_SieveRange; j += orderOfB)
         ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+j, -1);
      
      for (j = jHash[3]; j < ii_SieveRange; j += orderOfB)
         ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+j, -1);
         
      return;
   }
   
   // First giant step
   if (jHash[0] != HASH_NOT_FOUND)
      ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+jHash[0], +1);
      
   if (jHash[1] != HASH_NOT_FOUND)
      ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+jHash[1], +1);
   
   if (jHash[2] != HASH_NOT_FOUND)
      ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+jHash[2], -1);
   
   if (jHash[3] != HASH_NOT_FOUND)
      ip_CarolKyneaApp->ReportFactor(p, ii_SieveLow+jHash[3], -1);

   // Remaining giant steps
   // b <- 1/b^m (mod p)
   MpRes resB = mp.pow(mp.nToRes(inv_pb), ii_BabySteps);

   uint32_t nBase = ii_SieveLow;
   
   for (uint32_t step = 1; step < ii_GiantSteps; step++)
   {
      resA[0] = mp.mul(resA[0], resB);
      resA[1] = mp.mul(resA[1], resB);
      resA[2] = mp.mul(resA[2], resB);
      resA[3] = mp.mul(resA[3], resB);

      jHash[0] = ip_HashTable->Lookup(resA[0]);
      jHash[1] = ip_HashTable->Lookup(resA[1]);
      jHash[2] = ip_HashTable->Lookup(resA[2]);
      jHash[3] = ip_HashTable->Lookup(resA[3]);
   
      nBase += ii_BabySteps;
   
      if (jHash[0] != HASH_NOT_FOUND)
         ip_CarolKyneaApp->ReportFactor(p, nBase+jHash[0], +1);
         
      if (jHash[1] != HASH_NOT_FOUND)
         ip_CarolKyneaApp->ReportFactor(p, nBase+jHash[1], +1);
      
      if (jHash[2] != HASH_NOT_FOUND)
         ip_CarolKyneaApp->ReportFactor(p, nBase+jHash[2], -1);
      
      if (jHash[3] != HASH_NOT_FOUND)
         ip_CarolKyneaApp->ReportFactor(p, nBase+jHash[3], -1);
   }
  }

uint32_t  CarolKyneaWorker::BabySteps(uint64_t p, MpArith mp, MpRes resBase)
{
   MpRes    resBJ, resBJ0;
   uint32_t j;
   
   resBJ = resBJ0 = mp.pow(resBase, ii_MinN);
  
   for (j = 0; j < ii_BabySteps; j++)
   {
      ip_HashTable->Insert(resBJ, j);
      
      resBJ = mp.mul(resBJ, resBase);
            
      if (resBJ == resBJ0)
         return j+1;
   }

   return 0;
}

