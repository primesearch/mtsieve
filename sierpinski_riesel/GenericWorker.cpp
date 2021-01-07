/* GenericWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "GenericWorker.h"

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
   
   xfree(mBD);
   xfree(mCK);
   xfree(mBDCK);
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

   InitializeWorker();
  
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

void  GenericWorker::InitializeWorker(void)
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
   
   mCK = (MpResVec *) xmalloc(ii_SequenceCount*sizeof(MpResVec));
   mBDCK = (MpResVec *) xmalloc(ii_SubsequenceCount*sizeof(MpResVec));
   mBD = (MpResVec *) xmalloc(ii_BestQ*sizeof(MpResVec));

   for (idx=0; idx<4; idx++)
      ip_HashTable[idx] = new HashTable(ii_BabySteps);
}

void  GenericWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t lastPrime = 0;
#ifdef HAVE_GPU_WORKERS
   uint64_t minGpuPrime = ip_App->GetMinGpuPrime();
#endif
   uint64_t p[4];
   uint32_t b[4];

   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   if (il_SmallPrimeSieveLimit > 0 && il_SmallPrimeSieveLimit > *it)
   {
      lastPrime = ProcessSmallPrimes();

      if (il_SmallPrimeSieveLimit <= lastPrime)
      {
         ip_SierpinskiRieselApp->SetRebuildNeeded();
         
         // This will trigger a retest of some primes between il_SmallPrimeSieveLimit and lastPrime
         // using large primes logic.  This should only be a few primes, so no big performance hit.
         SetLargestPrimeTested(il_SmallPrimeSieveLimit, 0);
      }
      
      return;
   }

   b[0] = b[1] = b[2] = b[3] = ii_Base;

   while (it != iv_Primes.end())
   {
      p[0] = *it;
      it++;
      
      p[1] = *it;
      it++;
      
      p[2] = *it;
      it++;
      
      p[3] = *it;
      it++;
         
      DiscreteLogLargePrimes(b, p);
      
      SetLargestPrimeTested(p[3], 4);
      
      if (p[3] >= maxPrime)
         return;

#ifdef HAVE_GPU_WORKERS
      if (lastPrime < minGpuPrime && p[3] > minGpuPrime)
      {      
         ip_SierpinskiRieselApp->SetRebuildNeeded();
         return;
      }
#endif
   }
   
   if (il_GenericSeveLimit > 0 && il_GenericSeveLimit < p[3])
      ip_SierpinskiRieselApp->SetRebuildNeeded();
}

uint64_t  GenericWorker::ProcessSmallPrimes(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t p[4];
   uint32_t b[4];
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

      p[primeCount] = thePrime;
      b[primeCount] = (ii_Base % thePrime);
      
      primeCount++;
      
      if (primeCount == 4)
      {   
         DiscreteLogSmallPrimes(b, p);
         
         SetLargestPrimeTested(p[3], 4);
         
         primeCount = 0;
         
         if (p[3] >= maxPrime)
            return lastPrime;
      }
   }

   if (primeCount == 0)
      return lastPrime;
   
   actualCount = primeCount;
   
   while (primeCount < 4)
   {
      p[primeCount] = p[primeCount-1];
      b[primeCount] = b[primeCount-1];

      primeCount++;
   }
      
   DiscreteLogSmallPrimes(b, p);
   
   SetLargestPrimeTested(p[3], actualCount);
   
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
   MpResVec mb = mp.toMp(b);
   
   SetupDicreteLog(b, p, mp, mb);
   
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
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, 0, j), (verifiedFactors < 5));
               
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
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, 0, j), (verifiedFactors < 5));
               
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
            uint32_t nTerm = N_TERM(ssIdx, 0, j);
            
            verifiedFactors = 0;
            
            while (nTerm < ii_MaxN)
            {
               verifiedFactors++;
               
               // If the first few are valid, then assume that the rest are valid.  This will
               // speed up testing of small primes.
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_IDX(ssIdx), nTerm, (verifiedFactors < 5));
               
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
   MpResVec mb = mp.toMp(b);
   
   SetupDicreteLog(b, p, mp, mb);
      
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
               ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, 0, j), true);
               
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
            ip_SierpinskiRieselApp->ReportFactor(p[pIdx], SEQ_IDX(ssIdx), N_TERM(ssIdx, 0, j), true);
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
            ip_SierpinskiRieselApp->ReportFactor(p[0], SEQ_IDX(ssIdx), N_TERM(ssIdx, i, j), true);

         j = ip_HashTable[1]->Lookup(mBDCK[ssIdx][1]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[1], SEQ_IDX(ssIdx), N_TERM(ssIdx, i, j), true);

         j = ip_HashTable[2]->Lookup(mBDCK[ssIdx][2]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[2], SEQ_IDX(ssIdx), N_TERM(ssIdx, i, j), true);

         j = ip_HashTable[3]->Lookup(mBDCK[ssIdx][3]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p[3], SEQ_IDX(ssIdx), N_TERM(ssIdx, i, j), true);
      }
   }
}

void  GenericWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GenericWorker::TestMiniPrimeChunk not implemented");
}

// Compute a number of values that we need for the discrete log
void  GenericWorker::SetupDicreteLog(uint32_t *b, uint64_t *p, MpArithVec mp, MpResVec mb)
{
   uint32_t qIdx, seqIdx, ssIdx;
   uint64_t imod[4], umod[4], temp[4];

   imod[0] = invmod64(b[0], p[0]);
   imod[1] = invmod64(b[1], p[1]);
   imod[2] = invmod64(b[2], p[2]);
   imod[3] = invmod64(b[3], p[3]);

   MpResVec mI = mp.toMp(imod);
   mBM = mI;
   mBD[0] = mp.one();
     
   for (qIdx=1; qIdx<ii_BestQ; qIdx++)
   {
      mBD[qIdx] = mBM;
      
      mBM = mp.mul(mBM, mI);
   }
    
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      temp[0] = lmod64(-ip_Sequences[seqIdx].c, p[0]);
      temp[1] = lmod64(-ip_Sequences[seqIdx].c, p[1]);
      temp[2] = lmod64(-ip_Sequences[seqIdx].c, p[2]);
      temp[3] = lmod64(-ip_Sequences[seqIdx].c, p[3]);
      
      umod[0] = umod64(ip_Sequences[seqIdx].k, p[0]);
      umod[1] = umod64(ip_Sequences[seqIdx].k, p[1]);
      umod[2] = umod64(ip_Sequences[seqIdx].k, p[2]);
      umod[3] = umod64(ip_Sequences[seqIdx].k, p[3]);
      
      imod[0] = invmod64(umod[0], p[0]);
      imod[1] = invmod64(umod[1], p[1]);
      imod[2] = invmod64(umod[2], p[2]);
      imod[3] = invmod64(umod[3], p[3]);
      
      MpResVec mTemp = mp.toMp(temp);
      mI = mp.toMp(imod);

      mCK[seqIdx] = mp.mul(mTemp, mI);
   }

   // Compute -c/(k*b^d) (mod p) for each subsequence.
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      seqIdx = ip_Subsequences[ssIdx].seqIdx;
      qIdx = ip_Subsequences[ssIdx].q;

      mBDCK[ssIdx] = mp.mul(mBD[qIdx], mCK[seqIdx]);
   }
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
