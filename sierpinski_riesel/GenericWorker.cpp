/* GenericWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "GenericWorker.h"
#include "SierpinskiRieselApp.h"

#define SEQ_PTR(ssIdx)        (ip_Subsequences[(ssIdx)].seqPtr)
#define N_TERM(ssIdx, j)      ((ii_SieveLow+(j))*ii_BestQ + ip_Subsequences[(ssIdx)].q)

GenericWorker::GenericWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) theApp;
   ib_CanUseCIsOneLogic = srApp->CanUseCIsOneLogic();
   il_MaxK = srApp->GetMaxK();
   
   // Everything we need is done in the constuctor of the parent class
   ib_Initialized = true;
   
   ssHash = NULL;
   mBD = NULL;
   mBDCK = NULL;
}

void  GenericWorker::CleanUp(void)
{
   uint32_t idx;
   
   for (idx=0; idx<4; idx++)
      delete ip_HashTable[idx];
   
   xfree(ssHash);
   xfree(mBD);
   xfree(mBDCK);
}

void  GenericWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{
   uint64_t   max_k = 0;
   uint64_t   max_c = 0;
   uint64_t   c;
   seq_t     *seqPtr;

   ii_BestQ = bestQ;
   
   seqPtr = ip_FirstSequence;
   do
   {
      c = abs(seqPtr->c);
      
      max_k = MAX(max_k, seqPtr->k);
      max_c = MAX(max_c, c);

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   InitializeWorker();
  
   il_SmallPrimeSieveLimit = ip_SierpinskiRieselApp->GetSmallSievePrimeLimit();
}

void  GenericWorker::InitializeWorker(void)
{
   uint32_t idx;   
   uint32_t r = ii_MaxN/ii_BestQ - ii_MinN/ii_BestQ + 1;
   double babyStepFactor = 1.0; // DEFAULT_BABY_STEP_FACTOR from srsieve

   // In the worst case we will do do one table insertion and one mulmod for ii_BabySteps
   // baby steps, then s table lookups and s mulmods for ii_GiantSteps giant steps. The
   // average case depends on how many solutions are found and how early in the loop they
   // are found, which I don't know how to analyse. However for the worst case we just want
   // to minimise ii_BabySteps + s*ii_GiantSteps subject to ii_BabySteps*ii_GiantSteps >= r
   // which is when ii_BabySteps = sqrt(s*r).

   ii_GiantSteps = MAX(1, sqrt((double) r/ii_SubsequenceCount/babyStepFactor));
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
   
   mBDCK = (MpResVec *) xmalloc(ii_SubsequenceCount*sizeof(MpResVec));
   mBD = (MpResVec *) xmalloc(ii_BestQ*sizeof(MpResVec));

   for (idx=0; idx<4; idx++)
      ip_HashTable[idx] = new HashTable(ii_BabySteps);
}

void  GenericWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t lastPrime = 0;
   uint64_t ps[4];
   uint32_t bs[4];
   
   if (il_SmallPrimeSieveLimit > 0 && il_SmallPrimeSieveLimit > il_PrimeList[0])
   {
      lastPrime = ProcessSmallPrimes();

      // If we can switch to the "large" prime discrete log logic, then this is the time to do it.
      if (il_SmallPrimeSieveLimit <= lastPrime)
      {
         ip_SierpinskiRieselApp->SetRebuildNeeded();
         
         // This will trigger a retest of some primes between il_SmallPrimeSieveLimit and lastPrime
         // using large primes logic.  This should only be a few primes, so no big performance hit.
         SetLargestPrimeTested(il_SmallPrimeSieveLimit, 0);
      }
      
      return;
   }

#if defined(USE_OPENCL) || defined(USE_METAL)
   bool     switchToGPUWorkers = false;
   
   // If this is the scenario where a CPU worker was created (even though CpuWorkerCount = 0)
   // then switching to GPU workers will delete the CPU worker when it is no longer needed.
   if (ip_App->GetCpuWorkerCount() == 0 && ip_App->GetGpuWorkerCount() > 0)
   {
      if (maxPrime > ip_App->GetMinGpuPrime())
      {
         maxPrime = ip_App->GetMinGpuPrime();
         switchToGPUWorkers = true;
      }
   }
#endif

   bs[0] = bs[1] = bs[2] = bs[3] = ii_Base;

   for (uint32_t pIdx=0; pIdx<ii_WorkSize; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
         
      DiscreteLogLargePrimes(bs, ps);
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] >= maxPrime)
      {
#if defined(USE_OPENCL) || defined(USE_METAL)
         // This can only be true if we can switch to only running GPU workers.  Since this
         // is a CPU worker.  Its life effectively ends upon return from this method.
         if (switchToGPUWorkers)
         {
            ip_SierpinskiRieselApp->UseGpuWorkersUponRebuild();
            ip_SierpinskiRieselApp->SetRebuildNeeded();
         }
#endif
         
         return;
      }
   }

   // Determine if we can switch to the CisOne workers.  Note that if il_MaxK < il_MinGpuPrime
   // that the CPU workers will be used instead of the GPU workers
   if (ib_CanUseCIsOneLogic && ps[3] > il_MaxK && ii_SequenceCount == 1)
      ip_SierpinskiRieselApp->SetRebuildNeeded();
}

uint64_t  GenericWorker::ProcessSmallPrimes(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t ps[4];
   uint32_t bs[4];
   uint64_t thePrime, lastPrime = 0;
   uint32_t primeCount = 0;
   uint32_t actualCount = 0;
      
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx++)
   {
      lastPrime = thePrime = il_PrimeList[pIdx];
      
      if (ii_Base % thePrime  == 0)
         continue;

      ps[primeCount] = thePrime;
      bs[primeCount] = (ii_Base % thePrime);
      
      primeCount++;
      
      if (primeCount == 4)
      {   
         DiscreteLogSmallPrimes(bs, ps);
         
         SetLargestPrimeTested(ps[3], 4);
         
         primeCount = 0;
         
         if (ps[3] >= maxPrime)
            return lastPrime;
      }
   }

   if (primeCount == 0)
      return lastPrime;
   
   actualCount = primeCount;
   
   while (primeCount < 4)
   {
      ps[primeCount] = ps[primeCount-1];
      bs[primeCount] = bs[primeCount-1];

      primeCount++;
   }
      
   DiscreteLogSmallPrimes(bs, ps);
   
   SetLargestPrimeTested(ps[3], actualCount);
   
   return lastPrime;
}

void  GenericWorker::DiscreteLogSmallPrimes(uint32_t *b, uint64_t *p)
{
   uint32_t i, j, ssIdx;
   uint32_t pIdx;
   uint32_t orderOfB[4];
   uint32_t solutionCount;
   uint32_t firstSolution;
   uint32_t verifiedFactors;
   
   MpArithVec mp(p);
   MpResVec mb = mp.nToRes(b);
   
   SetupDiscreteLog(b, p, mp, mb);
   
   BabySteps(mp, mb, orderOfB);
   
   // b <- 1/b^m (mod p)
   if (ii_GiantSteps > 1)
      mBM = mp.pow(mBM, ii_BabySteps);

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
            j = ip_HashTable[pIdx]->Lookup(mBDCK[ssIdx][pIdx]);

            verifiedFactors = 0;
            
            while (j < ii_SieveRange)
            {
               verifiedFactors++;
               
               // If the first few are valid, then assume that the rest are valid.  This will
               // speed up testing of small primes.
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_PTR(ssIdx), N_TERM(ssIdx, j), (verifiedFactors < 5));
               
               j += orderOfB[pIdx];
            }
         }
         
         continue;
      }
      
      // First giant step
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         j = ip_HashTable[pIdx]->Lookup(mBDCK[ssIdx][pIdx]);
         
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
         for (i=1; i<ii_GiantSteps && solutionCount<=ii_SubsequenceCount; i++)
         {
            for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
            {
               if (ssHash[ssIdx] == HASH_NOT_FOUND || firstSolution == ssIdx)
               {
                  mBDCK[ssIdx][pIdx] = mp.mul(mBDCK[ssIdx], mBM, pIdx);

                  j = ip_HashTable[pIdx]->Lookup(mBDCK[ssIdx][pIdx]);
            
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
      }
   
      if (orderOfB[pIdx] > 0)
      {
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
         {
            j = ssHash[ssIdx];

            verifiedFactors = 0;
            
            while (j < ii_SieveRange)
            {
               verifiedFactors++;
               
               // If the first few are valid, then assume that the rest are valid.  This will
               // speed up testing of small primes.
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_PTR(ssIdx), N_TERM(ssIdx, j), (verifiedFactors < 5));
               
               j += orderOfB[pIdx];
            }
         }
      
         continue;
      }
      
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         j = ssHash[ssIdx];
         
         if (j != HASH_NOT_FOUND)
         {
            uint32_t nTerm = N_TERM(ssIdx, j);
            
            verifiedFactors = 0;
            
            while (nTerm <= ii_MaxN)
            {
               verifiedFactors++;
               
               // If the first few are valid, then assume that the rest are valid.  This will
               // speed up testing of small primes.
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_PTR(ssIdx), nTerm, (verifiedFactors < 5));
               
               nTerm += (p[pIdx] - 1);
            }
         }
      }
   }
}

void  GenericWorker::DiscreteLogLargePrimes(uint32_t *b, uint64_t *p)
{
   uint32_t i, j, ssIdx;
   uint32_t pIdx;
   uint32_t orderOfB[4];
   
   MpArithVec mp(p);
   MpResVec mb = mp.nToRes(b);
   
   SetupDiscreteLog(b, p, mp, mb);
      
   BabySteps(mp, mb, orderOfB);

   for (pIdx=0; pIdx<4; pIdx++)
   {
      if (orderOfB[pIdx] > 0)
      {
         for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
         {
            j = ip_HashTable[pIdx]->Lookup(mBDCK[ssIdx][pIdx]);
            
            while (j < ii_SieveRange)
            {
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_PTR(ssIdx), N_TERM(ssIdx, j), true);
               
               j += orderOfB[pIdx];
            }
         }
         
         continue;
      }

      // First giant step
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         j = ip_HashTable[pIdx]->Lookup(mBDCK[ssIdx][pIdx]);
         
         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_PTR(ssIdx), N_TERM(ssIdx, j), true);
      }
   }
   
   if (ii_GiantSteps < 2)
      return;

   // b <- 1/b^m (mod p)
   mBM = mp.pow(mBM, ii_BabySteps);
   
   for (i=1; i<ii_GiantSteps; i++)
   {
      for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
      {
         mBDCK[ssIdx] = mp.mul(mBDCK[ssIdx], mBM);

         j = ip_HashTable[0]->Lookup(mBDCK[ssIdx][0]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[0], SEQ_PTR(ssIdx), N_TERM(ssIdx, j + i*ii_BabySteps), true);

         j = ip_HashTable[1]->Lookup(mBDCK[ssIdx][1]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[1], SEQ_PTR(ssIdx), N_TERM(ssIdx, j + i*ii_BabySteps), true);

         j = ip_HashTable[2]->Lookup(mBDCK[ssIdx][2]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[2], SEQ_PTR(ssIdx), N_TERM(ssIdx, j + i*ii_BabySteps), true);

         j = ip_HashTable[3]->Lookup(mBDCK[ssIdx][3]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[3], SEQ_PTR(ssIdx), N_TERM(ssIdx, j + i*ii_BabySteps), true);
      }
   }
}

void  GenericWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GenericWorker::TestMiniPrimeChunk not implemented");
}

// Compute a number of values that we need for the discrete log
void  GenericWorker::SetupDiscreteLog(uint32_t *b, uint64_t *p, MpArithVec mp, MpResVec mb)
{
   uint32_t   qIdx, ssIdx;
   uint64_t   imod[4], umod[4], temp[4];
   seq_t     *seqPtr;

   imod[0] = invmod64(b[0], p[0]);
   imod[1] = invmod64(b[1], p[1]);
   imod[2] = invmod64(b[2], p[2]);
   imod[3] = invmod64(b[3], p[3]);

   MpResVec mCK;
   MpResVec mI = mp.nToRes(imod);
   mBM = mI;
   mBD[0] = mp.one();
     
   for (qIdx=1; qIdx<ii_BestQ; qIdx++)
   {
      mBD[qIdx] = mBM;
      
      mBM = mp.mul(mBM, mI);
   }
   
   seqPtr = ip_FirstSequence;
   do
   {
      temp[0] = lmod64(-seqPtr->c, p[0]);
      temp[1] = lmod64(-seqPtr->c, p[1]);
      temp[2] = lmod64(-seqPtr->c, p[2]);
      temp[3] = lmod64(-seqPtr->c, p[3]);
      
      umod[0] = umod64(seqPtr->k, p[0]);
      umod[1] = umod64(seqPtr->k, p[1]);
      umod[2] = umod64(seqPtr->k, p[2]);
      umod[3] = umod64(seqPtr->k, p[3]);
      
      imod[0] = invmod64(umod[0], p[0]);
      imod[1] = invmod64(umod[1], p[1]);
      imod[2] = invmod64(umod[2], p[2]);
      imod[3] = invmod64(umod[3], p[3]);
      
      MpResVec mTemp = mp.nToRes(temp);
      mI = mp.nToRes(imod);

      mCK = mp.mul(mTemp, mI);

      // Compute -c/(k*b^d) (mod p) for each subsequence.
      for (ssIdx=seqPtr->ssIdxFirst; ssIdx<=seqPtr->ssIdxLast; ssIdx++)
      {
         qIdx = ip_Subsequences[ssIdx].q;

         mBDCK[ssIdx] = mp.mul(mBD[qIdx], mCK);
      }

      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);

}

void  GenericWorker::BabySteps(MpArithVec mp, MpResVec mb, uint32_t *orderOfB)
{
   uint32_t j, pIdx;
   uint64_t resBJ;
   
   const MpResVec mBexpQ = mp.pow(mb, ii_BestQ);
   const MpResVec mBJ = mp.pow(mBexpQ, ii_SieveLow);
      
   for (pIdx=0; pIdx<4; pIdx++)
   {
      ip_HashTable[pIdx]->Clear();
      
      orderOfB[pIdx] = 0;
      resBJ = mBJ[pIdx];
            
      for (j=0; j<ii_BabySteps; j++)
      {
         ip_HashTable[pIdx]->Insert(resBJ, j);
         
         resBJ = mp.mul(resBJ, mBexpQ, pIdx);
         
         if (resBJ == mBJ[pIdx])
         {
            orderOfB[pIdx] = j + 1;
            break;
         }
      }
   }
}
