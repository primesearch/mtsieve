/* LegendreWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm/sse-asm-x86.h"
#include "LegendreWorker.h"

#define SEQ_IDX(ssIdx)     (ip_WorkerSubsequences[(ssIdx)].seqIdx)
#define N_TERM(ssIdx, m)   ((m)*ii_BestQ + ip_WorkerSubsequences[ssIdx].d)

LegendreWorker::LegendreWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  LegendreWorker::CleanUp(void)
{
   uint32_t idx;

   xfree(ssHash);
   
   for (idx=0; idx<4; idx++)
      delete ip_HashTable[idx];
   
   for (idx=0; idx<ii_BestQ; idx++)
      xfree(bd64[idx]);   

   xfree(bd64);
   
   for (idx=0; idx<ii_SubsequenceCount; idx++)
      xfree(bdck64[idx]);   

   xfree(bdck64);
   
   for (idx=0; idx<ii_SequenceCount; idx++)
      xfree(ck64[idx]);   

   xfree(ck64);
   
   xfree(ip_WorkerSubsequences);
   xfree(ip_WorkerSequences);
}

void  LegendreWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t primeList[4];
   
   vector<uint64_t>::iterator it = iv_Primes.begin();
   
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

      DiscreteLog(primeList);
      
      SetLargestPrimeTested(primeList[3], 4);
      
      if (primeList[3] >= maxPrime)
         break;
   }
}

void  LegendreWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("LegendreWorker::TestMiniPrimeChunk not implemented");
}

void  LegendreWorker::SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, helper_seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount)
{
   if (ip_WorkerSequences)
      CleanUp();
   
   ii_BestQ = bestQ;
   ii_SequenceCount = sequenceCount;
   ip_WorkerSequences = (worker_seq_t *) xmalloc(sequenceCount * sizeof(worker_seq_t));
   
   for (uint32_t i=0; i<sequenceCount; i++)
   {
      ip_WorkerSequences[i].k = sequences[i].k;
      ip_WorkerSequences[i].c = sequences[i].c;
   }
   
   ii_SubsequenceCount = subsequenceCount;
   ip_WorkerSubsequences = (worker_subseq_t *) xmalloc(subsequenceCount * sizeof(worker_subseq_t));
   
   for (uint32_t i=0; i<subsequenceCount; i++)
   {
      ip_WorkerSubsequences[i].seqIdx = subsequences[i].seqIdx;
      ip_WorkerSubsequences[i].k = subsequences[i].k;
      ip_WorkerSubsequences[i].c = subsequences[i].c;
      ip_WorkerSubsequences[i].d = subsequences[i].d;
   }
   
   InitializeDiscreteLog();
}

void  LegendreWorker::InitializeDiscreteLog(void)
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

void  LegendreWorker::DiscreteLog(uint64_t *primeList)
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
   
   BuildTablesSSE(b, primeList, bm64);
   
   // b <- base^Q (mod p)
   sse_powmod_4b_1n_4p(b, ii_BestQ, primeList, invp);
   
   bj0[0] = b[0];
   bj0[1] = b[1];
   bj0[2] = b[2];
   bj0[3] = b[3];
   
   sse_powmod_4b_1n_4p(bj0, ii_SieveLow, primeList, invp);
   
   BabyStepsSSE(primeList, invp, b, bj0, orderOfB);

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

void  LegendreWorker::BabyStepsSSE(uint64_t *primeList, double *invp, uint64_t *b, uint64_t *bj0, uint32_t *orderOfB)
{
   uint32_t j, pIdx;
   uint64_t bj[4];
   uint64_t bs[4];
   uint32_t have = 0;
   bool     flag[4];
   double   bdivp[4];
   
   for (pIdx=0; pIdx<4; pIdx++)
   {
      ip_HashTable[pIdx]->Clear();
      
      orderOfB[pIdx] = 0;
      bj[pIdx] = bj0[pIdx];
      bs[pIdx] = b[pIdx];
      flag[pIdx] = false;
      bdivp[pIdx] = invp[pIdx] * (double) bs[pIdx];
   }

   for (j=0; j<ii_BabySteps; j++)
   {
      if (!flag[0]) ip_HashTable[0]->Insert(bj[0], j);
      if (!flag[1]) ip_HashTable[1]->Insert(bj[1], j);
      if (!flag[2]) ip_HashTable[2]->Insert(bj[2], j);
      if (!flag[3]) ip_HashTable[3]->Insert(bj[3], j);
    
      sse_mulmod_4a_4b_4p(bj, bs, primeList, bdivp);

      for (pIdx=0; pIdx<4; pIdx++)
         if (!flag[pIdx] && bj[pIdx] == bj0[pIdx])
         {
            orderOfB[pIdx] = j + 1;
            flag[pIdx] = true;
            
            if (++have == 4)
               return;
         }
   }
}

void   LegendreWorker::BuildTablesSSE(uint64_t *baseToUse, uint64_t *primeList, uint64_t *bm64)
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
      ck64[seqIdx][0] = smod64(-ip_WorkerSequences[seqIdx].c, primeList[0]);
      ck64[seqIdx][1] = smod64(-ip_WorkerSequences[seqIdx].c, primeList[1]);
      ck64[seqIdx][2] = smod64(-ip_WorkerSequences[seqIdx].c, primeList[2]);
      ck64[seqIdx][3] = smod64(-ip_WorkerSequences[seqIdx].c, primeList[3]);
      
      umod[0] = umod64(ip_WorkerSequences[seqIdx].k, primeList[0]);
      umod[1] = umod64(ip_WorkerSequences[seqIdx].k, primeList[1]);
      umod[2] = umod64(ip_WorkerSequences[seqIdx].k, primeList[2]);
      umod[3] = umod64(ip_WorkerSequences[seqIdx].k, primeList[3]);
      
      inv[0] = invmod64(umod[0], primeList[0]);
      inv[1] = invmod64(umod[1], primeList[1]);
      inv[2] = invmod64(umod[2], primeList[2]);
      inv[3] = invmod64(umod[3], primeList[3]);
      
      fpu_mulmod_4a_4b_4p(ck64[seqIdx], inv, primeList);
   }

   // Compute -c/(k*b^d) (mod p) for each subsequence.
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      bdck64[ssIdx][0] = ck64[ip_WorkerSubsequences[ssIdx].seqIdx][0];
      bdck64[ssIdx][1] = ck64[ip_WorkerSubsequences[ssIdx].seqIdx][1];
      bdck64[ssIdx][2] = ck64[ip_WorkerSubsequences[ssIdx].seqIdx][2];
      bdck64[ssIdx][3] = ck64[ip_WorkerSubsequences[ssIdx].seqIdx][3];
      
      fpu_mulmod_4a_4b_4p(bdck64[ssIdx], bd64[ip_WorkerSubsequences[ssIdx].d], primeList);
   }
   
   fpu_pop();
   fpu_pop();
   fpu_pop();
   fpu_pop();
}
