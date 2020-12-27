/* GenericWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "GenericWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"

GenericWorker::GenericWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   // Everything we need is done in the constuctor of the parent class
   ib_Initialized = true;
}

void  GenericWorker::CleanUp(void)
{
   uint32_t idx;

   xfree(ssHash);
   
   for (idx=0; idx<4; idx++)
      delete ip_HashTable[idx];
   
   for (idx=0; idx<ii_BestQ; idx++)
   {
      xfree(bd64[idx]);
   }

   xfree(bd64);
   
   for (idx=0; idx<ii_SubsequenceCount; idx++)
   {
      xfree(bdck64[idx]);
   }

   xfree(bdck64);
   
   for (idx=0; idx<ii_SequenceCount; idx++)
   {
      xfree(ck64[idx]);
   }

   xfree(ck64);
}

void  GenericWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t primeList[4];

   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   if (il_SmallPrimeSieveLimit > 0 && il_SmallPrimeSieveLimit > *it)
   {
      maxPrime = ProcessSmallPrimes();

      if (il_SmallPrimeSieveLimit <= maxPrime)
      {
         ip_SierpinskiRieselApp->SetRebuildNeeded();
         
         // This will trigger a retest of some primes between il_SmallPrimeSieveLimit and maxPrime
         // using large primes logic.  This should only be a few primes, so no big performance hit.
         SetLargestPrimeTested(il_SmallPrimeSieveLimit, 0);
      }
      
      return;
   }
    
   while (it != iv_Primes.end())
   {
      primeList[0] = *it;
      it++;
      
      primeList[1] = *it;
      it++;
      
      primeList[2] = *it;
      it++;
      
      primeList[3] = *it;
      it++;

      DiscreteLogLargePrimes(primeList);
      
      SetLargestPrimeTested(primeList[3], 4);
      
      if (primeList[3] >= maxPrime)
         return;
   }
   
   if (il_GenericSeveLimit > 0 && il_GenericSeveLimit < primeList[3])
      ip_SierpinskiRieselApp->SetRebuildNeeded();
}

void  GenericWorker::SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount)
{
   uint64_t max_k = 0;
   uint64_t max_c = 0;
   uint64_t c;

   if (ip_Sequences)
      CleanUp();

   ii_BestQ = bestQ;

   ip_Sequences = sequences;
   ii_SequenceCount = sequenceCount;
   
   for (uint32_t seqIdx=0; seqIdx<sequenceCount; seqIdx++)
   {
      c = abs(ip_Sequences[seqIdx].c);
      
      max_k = MAX(max_k, ip_Sequences[seqIdx].k);
      max_c = MAX(max_c, c);
   }
   
   ip_Subsequences = subsequences;
   ii_SubsequenceCount = subsequenceCount;

   InitializeDiscreteLog();
  
   il_SmallPrimeSieveLimit = ip_SierpinskiRieselApp->GetSmallSievePrimeLimit();
   
   if (largestPrimeTested > il_SmallPrimeSieveLimit)
   {
      il_SmallPrimeSieveLimit = 0;

      // We also have the possibility of switching to the CIsOne worker, but can only
      // do so if max(c) = 1 and we have sieved past max(k) and max(n).
      if (largestPrimeTested < MAX(max_k, ii_MaxN) && max_c == 1)
         il_GenericSeveLimit = MAX(max_k, ii_MaxN);
      else
         il_GenericSeveLimit = 0;
   }
   
   // The AVX code doesn't work and isn't much faster, so it isn't worth using.
   //if (ip_SierpinskiRieselApp->UseAvxIfAvailable() && CpuSupportsAvx())
   //   SetMiniChunkRange(il_SmallPrimeSieveLimit+1, PMAX_MAX_52BIT, AVX_ARRAY_SIZE);
}

void  GenericWorker::InitializeDiscreteLog(void)
{
   uint32_t idx;   
   uint32_t r = ii_MaxN/ii_BestQ - ii_MinN/ii_BestQ + 1;
   double baby_step_factor = 1.0; // DEFAULT_BABY_STEP_FACTOR from srsieve

   // In the worst case we will do do one table insertion and one mulmod for ii_BabySteps
   // baby steps, then s table lookups and s mulmods for ii_GiantSteps giant steps. The
   // average case depends on how many solutions are found and how early in the loop they
   // are found, which I don't know how to analyse. However for the worst case we just want
   // to minimise ii_BabySteps + s*ii_GiantSteps subject to ii_BabySteps*ii_GiantSteps >= r
   // which is when ii_BabySteps = sqrt(s*r).

   ii_GiantSteps = MAX(1, sqrt((double) r/ii_SubsequenceCount/baby_step_factor));
   ii_BabySteps = MIN(r, ceil((double) r/ii_GiantSteps));

   if (ii_BabySteps > HASH_MAX_ELTS)
   {
      ii_GiantSteps = ceil((double)r/HASH_MAX_ELTS);
      ii_BabySteps = ceil((double)r/ii_GiantSteps);
   }

   ii_SieveLow = ii_MinN / ii_BestQ;
   ii_SieveRange = ii_BabySteps*ii_GiantSteps;

   if (ii_SieveLow > ii_MinN/ii_BestQ)
      FatalError("ii_SieveLow was not computed correctly");
   
   if (ii_MaxN/ii_BestQ >= ii_SieveLow+ii_SieveRange)
      FatalError("ii_SieveRange was not computed correctly");
   ssHash = (uint32_t *) xmalloc(ii_SubsequenceCount*sizeof(uint32_t *));
   
   ck64 = (uint64_t **) xmalloc(ii_SequenceCount*sizeof(uint64_t *));
   
   for (idx=0; idx<ii_SequenceCount; idx++)
   {
      ck64[idx] = (uint64_t *) xmalloc(4*sizeof(uint64_t));
   }
   
   bdck64 = (uint64_t **) xmalloc(ii_SubsequenceCount*sizeof(uint64_t *));
   
   for (idx=0; idx<ii_SubsequenceCount; idx++)
   {
      bdck64[idx] = (uint64_t *) xmalloc(4*sizeof(uint64_t));
   }
   
   bd64 = (uint64_t **) xmalloc(ii_BestQ*sizeof(uint64_t *));
   
   for (idx=0; idx<ii_BestQ; idx++)
   {
      bd64[idx] = (uint64_t *) xmalloc(4*sizeof(uint64_t));
   }
   
   for (idx=0; idx<4; idx++)
      ip_HashTable[idx] = new HashTable(ii_BabySteps);
}

uint64_t  GenericWorker::ProcessSmallPrimes(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t primeList[4];
   uint64_t bases[4];
   uint64_t thePrime, lastPrime = 0;
   uint32_t primeCount = 0;
   uint32_t actualCount = 0;
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   while (it != iv_Primes.end())
   {
      lastPrime = thePrime = *it;
      it++;
      
      if (ii_Base % thePrime  == 0)
         continue;

      primeList[primeCount] = thePrime;
      bases[primeCount] = (ii_Base % thePrime);
      
      primeCount++;
      
      if (primeCount == 4)
      {
         DiscreteLogSmallPrimes(primeList, bases);
         
         SetLargestPrimeTested(primeList[3], 4);
         
         primeCount = 0;
         
         if (primeList[3] >= maxPrime)
            return lastPrime;
      }
   }

   if (primeCount == 0)
      return lastPrime;
   
   actualCount = primeCount;
   
   while (primeCount < 4)
   {
      primeList[primeCount] = primeList[primeCount-1];
      bases[primeCount] = bases[primeCount-1];

      primeCount++;
   }
   
   DiscreteLogSmallPrimes(primeList, bases);
   
   SetLargestPrimeTested(primeList[3], actualCount);
   
   return lastPrime;
}

void  GenericWorker::DiscreteLogSmallPrimes(uint64_t *primeList, uint64_t *bases)
{
   uint32_t i, j, ssIdx;
   uint32_t pIdx;
   uint32_t orderOfB[4];
   uint32_t solutionCount;
   uint32_t firstSolution;
   uint64_t baseList[4], b[4], bm64[4], bj0[4];
   double   invp[4];
   
   b[0] = baseList[0] = bases[0];
   b[1] = baseList[1] = bases[1];
   b[2] = baseList[2] = bases[2];
   b[3] = baseList[3] = bases[3];
  
   invp[0] = 1.0 / (double) primeList[0];
   invp[1] = 1.0 / (double) primeList[1];
   invp[2] = 1.0 / (double) primeList[2];
   invp[3] = 1.0 / (double) primeList[3];
      
   BuildTables(baseList, primeList, invp, bm64);
   
   // b <- base^Q (mod p)
   sse_powmod_4b_1n_4p(b, ii_BestQ, primeList, invp);
     
   bj0[0] = b[0];
   bj0[1] = b[1];
   bj0[2] = b[2];
   bj0[3] = b[3];
   
   sse_powmod_4b_1n_4p(bj0, ii_SieveLow, primeList, invp);

   BabyStepsSmallPrimes(primeList, b, bj0, orderOfB);
   
   if (ii_GiantSteps > 1)
      // b <- 1/b^m (mod p)
      sse_powmod_4b_1n_4p(bm64, ii_BabySteps, primeList, invp);
      //fpu_powmod_4b_1n_4p(bm64, ii_BabySteps, primeList);
   
   for (pIdx=0; pIdx<4; pIdx++)
   {      
      solutionCount = 0;
      firstSolution = HASH_NOT_FOUND;
      
      for (i=0; i<ii_SubsequenceCount; i++)
         ssHash[i] = HASH_NOT_FOUND;
            
      if (orderOfB[pIdx] > 0)
      {         
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
         {
            j = ip_HashTable[pIdx]->Lookup(bdck64[ssIdx][pIdx]);

            while (j < ii_SieveRange)
            {
               ip_SierpinskiRieselApp->ReportFactor(primeList[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, ii_SieveLow+j));
               
               j += orderOfB[pIdx];
            }
         }
         
         continue;
      }
      
      // First giant step
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         j = ip_HashTable[pIdx]->Lookup(bdck64[ssIdx][pIdx]);
         
         if (j != HASH_NOT_FOUND)
         {
            solutionCount++;
            ssHash[ssIdx] = j;
            
            if (solutionCount == 1)
               firstSolution = ssIdx;
         }
      }

      // Remaining giant steps
      if (ii_GiantSteps > 1)
      {
         fpu_push_adivb(bm64[pIdx], primeList[pIdx]);
         
         for (i=1; i<ii_GiantSteps && solutionCount<=ii_SubsequenceCount; i++)
         {
            for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
            {
               if (ssHash[ssIdx] == HASH_NOT_FOUND || firstSolution == ssIdx)
               {               
                  bdck64[ssIdx][pIdx] = fpu_mulmod_iter(bdck64[ssIdx][pIdx], bm64[pIdx], primeList[pIdx]);

                  j = ip_HashTable[pIdx]->Lookup(bdck64[ssIdx][pIdx]);
            
                  if (j != HASH_NOT_FOUND)
                  {
                     solutionCount++;
                     
                     if (firstSolution == ssIdx) /* repeat solution */
                     {
                        orderOfB[pIdx] = i*ii_BabySteps+j - ssHash[ssIdx];
                        firstSolution = HASH_NOT_FOUND;   /* no more repeats needed */
                     }
                     else
                     {
                        ssHash[ssIdx] = i*ii_BabySteps+j;
                        
                        if (solutionCount == 1) /* first solution */
                           firstSolution = ssIdx;
                     }
                  }
               }
            }
         }

         fpu_pop();
      }
   
      if (orderOfB[pIdx] > 0)
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
            for (j = ssHash[ssIdx]; j < ii_SieveRange; j += orderOfB[pIdx])
               ip_SierpinskiRieselApp->ReportFactor(primeList[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, ii_SieveLow+j));
      else
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
            if (ssHash[ssIdx] != HASH_NOT_FOUND)
            {
               uint32_t nTerm = N_TERM(ssIdx, ii_SieveLow+ssHash[ssIdx]);
               
               while (nTerm < ii_MaxN)
               {   
                  ip_SierpinskiRieselApp->ReportFactor(primeList[pIdx], SEQ_IDX(ssIdx), nTerm);
                  nTerm += (primeList[pIdx] - 1);
               }
            }
   }
}

void  GenericWorker::DiscreteLogLargePrimes(uint64_t *primeList)
{
   uint32_t i, j, ssIdx;
   uint32_t pIdx;
   uint32_t orderOfB[4];   
   uint64_t b[4], bm64[4], bj0[4];
   double   invp[4];
      
   b[0] = ii_Base;
   b[1] = ii_Base;
   b[2] = ii_Base;
   b[3] = ii_Base;
   
   invp[0] = 1.0 / (double) primeList[0];
   invp[1] = 1.0 / (double) primeList[1];
   invp[2] = 1.0 / (double) primeList[2];
   invp[3] = 1.0 / (double) primeList[3];
   
   BuildTables(b, primeList, invp, bm64);
   
   // b <- base^Q (mod p)
   sse_powmod_4b_1n_4p(b, ii_BestQ, primeList, invp);
   
   bj0[0] = b[0];
   bj0[1] = b[1];
   bj0[2] = b[2];
   bj0[3] = b[3];
   
   // This powmod routine doesn't support n = 0 and it is possible for ii_SieveLow to be 0
   if (ii_SieveLow == 0)
      bj0[0] = bj0[1] = bj0[2] = bj0[3] = 1;
   else
      sse_powmod_4b_1n_4p(bj0, ii_SieveLow, primeList, invp);
   
   BabyStepsBigPrimes(primeList, invp, b, bj0, orderOfB);

   for (pIdx=0; pIdx<4; pIdx++)
   {
      if (orderOfB[pIdx] > 0)
      {
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
         {
            j = ip_HashTable[pIdx]->Lookup(bdck64[ssIdx][pIdx]);
            
            while (j < ii_SieveRange)
            {
               ip_SierpinskiRieselApp->ReportFactor(primeList[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, ii_SieveLow+j));
               
               j += orderOfB[pIdx];
            }
         }
         
         continue;
      }

      // First giant step
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         j = ip_HashTable[pIdx]->Lookup(bdck64[ssIdx][pIdx]);
         
         if (j != HASH_NOT_FOUND)               
            ip_SierpinskiRieselApp->ReportFactor(primeList[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, ii_SieveLow+j));
      }
   }
   
   if (ii_GiantSteps < 2)
      return;
   
   // b <- 1/b^m (mod p)
   sse_powmod_4b_1n_4p(bm64, ii_BabySteps, primeList, invp);
   //fpu_powmod_4b_1n_4p(bm64, ii_BabySteps, primeList);
   
   fpu_push_1divp(primeList[3]);
   fpu_push_1divp(primeList[2]);
   fpu_push_1divp(primeList[1]);
   fpu_push_1divp(primeList[0]);
      
   for (i=1; i<ii_GiantSteps; i++)
   {
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         fpu_mulmod_4a_4b_4p(bdck64[ssIdx], bm64, primeList);

         for (pIdx=0; pIdx<4; pIdx++)
         {
            j = ip_HashTable[pIdx]->Lookup(bdck64[ssIdx][pIdx]);

            if (j != HASH_NOT_FOUND)
               ip_SierpinskiRieselApp->ReportFactor(primeList[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, ii_SieveLow+i*ii_BabySteps+j));
         }
      }
   }
      
   fpu_pop();
   fpu_pop();
   fpu_pop();
   fpu_pop();
}

void  GenericWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GenericWorker::TestMiniPrimeChunk not implemented");
}

void  GenericWorker::BabyStepsSmallPrimes(uint64_t *primeList, uint64_t *b, uint64_t *bj0, uint32_t *orderOfB)
{
   uint32_t j, pIdx;
   uint64_t bj;
   
   for (pIdx=0; pIdx<4; pIdx++)
   {
      ip_HashTable[pIdx]->Clear();
      
      orderOfB[pIdx] = 0;
      bj = bj0[pIdx];
      
      fpu_push_adivb(b[pIdx], primeList[pIdx]);
      
      for (j=0; j<ii_BabySteps; j++)
      {         
         ip_HashTable[pIdx]->Insert(bj, j);
         
         bj = fpu_mulmod_iter(bj, b[pIdx], primeList[pIdx]);
         
         if (bj == bj0[pIdx])
         {
            orderOfB[pIdx] = j + 1;
            break;
         }
      }
      
      fpu_pop();
   }
}

void  GenericWorker::BabyStepsBigPrimes(uint64_t *primeList, double *invp, uint64_t *b, uint64_t *bj0, uint32_t *orderOfB)
{
   uint32_t j, pIdx;
   uint64_t bj;
   double   binvp;
  
   for (pIdx=0; pIdx<4; pIdx++)
   {
      ip_HashTable[pIdx]->Clear();
      
      orderOfB[pIdx] = 0;
      bj = bj0[pIdx];
      
      binvp = invp[pIdx] * (double) b[pIdx];
            
      for (j=0; j<ii_BabySteps; j++)
      {         
         ip_HashTable[pIdx]->Insert(bj, j);
         
         bj = sse_mulmod(bj, b[pIdx], primeList[pIdx], &binvp);
         
         if (bj == bj0[pIdx])
         {
            orderOfB[pIdx] = j + 1;
            break;
         }
      }
   }
}

void   GenericWorker::BuildTables(uint64_t *baseToUse, uint64_t *primeList, double *invp, uint64_t *bm64)
{
   uint64_t inv_b[4];
   uint32_t qIdx, seqIdx, ssIdx;
   uint64_t umod[4], inv[4];

   fpu_push_1divp(primeList[3]);
   fpu_push_1divp(primeList[2]);
   fpu_push_1divp(primeList[1]);
   fpu_push_1divp(primeList[0]);
   
   // Precompute 1/b^d (mod p) for 0 <= d <= Q.

   bd64[0][0] = bd64[0][1] = bd64[0][2] = bd64[0][3] = 1;
   
   bm64[0] = inv_b[0] = invmod64(baseToUse[0], primeList[0]);
   bm64[1] = inv_b[1] = invmod64(baseToUse[1], primeList[1]);
   bm64[2] = inv_b[2] = invmod64(baseToUse[2], primeList[2]);
   bm64[3] = inv_b[3] = invmod64(baseToUse[3], primeList[3]);

   for (qIdx=1; qIdx<ii_BestQ; qIdx++)
   {
      bd64[qIdx][0] = bm64[0];
      bd64[qIdx][1] = bm64[1];
      bd64[qIdx][2] = bm64[2];
      bd64[qIdx][3] = bm64[3];
         
      fpu_mulmod_4a_4b_4p(bm64, inv_b, primeList);
   }
    
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      ck64[seqIdx][0] = smod64(-ip_Sequences[seqIdx].c, primeList[0]);
      ck64[seqIdx][1] = smod64(-ip_Sequences[seqIdx].c, primeList[1]);
      ck64[seqIdx][2] = smod64(-ip_Sequences[seqIdx].c, primeList[2]);
      ck64[seqIdx][3] = smod64(-ip_Sequences[seqIdx].c, primeList[3]);
      
      umod[0] = umod64(ip_Sequences[seqIdx].k, primeList[0]);
      umod[1] = umod64(ip_Sequences[seqIdx].k, primeList[1]);
      umod[2] = umod64(ip_Sequences[seqIdx].k, primeList[2]);
      umod[3] = umod64(ip_Sequences[seqIdx].k, primeList[3]);
      
      inv[0] = invmod64(umod[0], primeList[0]);
      inv[1] = invmod64(umod[1], primeList[1]);
      inv[2] = invmod64(umod[2], primeList[2]);
      inv[3] = invmod64(umod[3], primeList[3]);
      
      fpu_mulmod_4a_4b_4p(ck64[seqIdx], inv, primeList);
   }

   // Compute -c/(k*b^d) (mod p) for each subsequence.
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      bdck64[ssIdx][0] = ck64[ip_Subsequences[ssIdx].seqIdx][0];
      bdck64[ssIdx][1] = ck64[ip_Subsequences[ssIdx].seqIdx][1];
      bdck64[ssIdx][2] = ck64[ip_Subsequences[ssIdx].seqIdx][2];
      bdck64[ssIdx][3] = ck64[ip_Subsequences[ssIdx].seqIdx][3];
      
      fpu_mulmod_4a_4b_4p(bdck64[ssIdx], bd64[ip_Subsequences[ssIdx].q], primeList);
   }
   
   fpu_pop();
   fpu_pop();
   fpu_pop();
   fpu_pop();
}