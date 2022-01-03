/* CisOneWithMultipleSequencesWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "CisOneWithMultipleSequencesWorker.h"
#include "../core/inline.h"
#include "../core/MpArith.h"

#define N_TERM(q, i, j)      ((ii_SieveLow + (j) + (i)*babySteps)*ii_BestQ + q)

CisOneWithMultipleSequencesWorker::CisOneWithMultipleSequencesWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ip_CisOneHelper = (CisOneWithMultipleSequencesHelper *) appHelper;
   
   ii_BaseMultiple = ip_CisOneHelper->GetBaseMultiple();
   ii_LimitBase = ip_CisOneHelper->GetLimitBase();
   ii_PowerResidueLcm = ip_CisOneHelper->GetPowerResidueLcm();
   
   // Everything we need is done in the constuctor of the parent class
   ib_Initialized = true;
}

void  CisOneWithMultipleSequencesWorker::CleanUp(void)
{
   delete ip_HashTable;
   
   xfree(ip_UsableSubsequences);
   
   xfree(resBD);
   xfree(resBJ);
   xfree(resX);
}

void  CisOneWithMultipleSequencesWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{
   ii_BestQ = bestQ;   
   ii_SieveLow = ii_MinN / ii_BestQ;

   resX = (MpRes *) xmalloc((ii_PowerResidueLcm+5) * sizeof(MpRes));
   resBD = (MpRes *) xmalloc((ii_SubsequenceCount+4)*sizeof(MpRes));
   resBJ = (MpRes *) xmalloc((ii_SubsequenceCount+4)*sizeof(MpRes));
   
   ip_UsableSubsequences = (useable_subseq_t *) xmalloc(ii_SubsequenceCount*sizeof(useable_subseq_t));
   
   ip_DivisorShifts = ip_CisOneHelper->GetDivisorShifts();
   ip_PowerResidueIndices = ip_CisOneHelper->GetPowerResidueIndices();
   ii_MaxBabySteps = ip_CisOneHelper->GetMaxBabySteps();
   
   ii_Dim1 = ip_CisOneHelper->GetDim1();
   ii_Dim2 = ip_CisOneHelper->GetDim2();
   ii_Dim3 = ip_CisOneHelper->GetDim3();
   
   ip_CongruentSubseqIndices = ip_CisOneHelper->GetCongruentSubseqIndices();
   ip_AllSubseqs = ip_CisOneHelper->GetAllSubseqs();
   ip_AllLadders = ip_CisOneHelper->GetAllLadders();   
   
   ip_Legendre = ip_CisOneHelper->GetLegendre();
   ip_LegendreTable = ip_CisOneHelper->GetLegendreTable();
   
   ip_HashTable = new HashTable(ii_MaxBabySteps);

   ib_AllSequencesHaveLegendreTables = true;
   for (uint32_t seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      legendre_t *legendrePtr = &ip_Legendre[seqIdx];
      
      if (!legendrePtr->haveMap)
         ib_AllSequencesHaveLegendreTables = false;
   }
}

void  CisOneWithMultipleSequencesWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t p;

#ifdef HAVE_GPU_WORKERS
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

   vector<uint64_t>::iterator it = iv_Primes.begin();
   
   while (it != iv_Primes.end())
   {
      p = *it;
      it++;

      TestSinglePrime(p);

      SetLargestPrimeTested(p, 1);
      
      if (p >= maxPrime)
      {
#ifdef HAVE_GPU_WORKERS
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
}
void  CisOneWithMultipleSequencesWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CisOneWithMultipleSequencesWorker::TestMiniPrimeChunk not implemented");
}

void  CisOneWithMultipleSequencesWorker::TestSinglePrime(uint64_t p)
{
   uint64_t invBase;
   uint32_t k, orderOfB, ssCount;
   uint32_t babySteps, giantSteps;
   uint32_t i, j;
   
   MpArith mp(p);

   // compute 1/base (mod p)
   invBase = invmod64(ii_Base, p);
      
   MpRes resBase = mp.nToRes(ii_Base);
   MpRes resInvBase = mp.nToRes(invBase);

   // -ckb^d is an r-th power residue for at least one term (k*b^d)*(b^Q)^(n/Q)+c of this subsequence
   ClimbLadder(mp, resBase);

   ssCount = SetupDiscreteLog(mp, resBase, resInvBase);

   // If no subsequences for this p, then no factors, so return
   if (ssCount == 0)
      return;
   
   ip_HashTable->Clear();
   
   babySteps = ip_Subsequences[ssCount-1].babySteps;
   giantSteps = ip_Subsequences[ssCount-1].giantSteps;

   orderOfB = BabySteps(mp, resBase, resInvBase, ii_MaxBabySteps);
     
   if (orderOfB > 0)
   {
      // If orderOfB > 0, then this is all the information we need to
      // determine every solution for this p, so no giant steps are neede
      for (k=0; k<ssCount; k++)
      {
          j = ip_HashTable->Lookup(ip_UsableSubsequences[k].resBDCK);

          while (j < babySteps * giantSteps)
          {
             ip_SierpinskiRieselApp->ReportFactor(p, ip_UsableSubsequences[k].seqPtr, N_TERM(ip_UsableSubsequences[k].q, 0, j), true);
             
             j += orderOfB;
          }
      }
   }
   else
   {
      // First giant step
      for (k=0; k<ssCount; k++)
      {
         j = ip_HashTable->Lookup(ip_UsableSubsequences[k].resBDCK);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p, ip_UsableSubsequences[k].seqPtr, N_TERM(ip_UsableSubsequences[k].q, 0, j), true);
      }

      // Remaining giant steps
      if (giantSteps > 1)
      {
         MpRes resBQM = mp.pow(resBexpQ, babySteps);
  
         for (i=1; i<giantSteps; i++)
         {
            for (k=0; k<ssCount; k++)
            {         
               ip_UsableSubsequences[k].resBDCK = mp.mul(ip_UsableSubsequences[k].resBDCK, resBQM);
               
               j = ip_HashTable->Lookup(ip_UsableSubsequences[k].resBDCK);

               if (j != HASH_NOT_FOUND)
                  ip_SierpinskiRieselApp->ReportFactor(p, ip_UsableSubsequences[k].seqPtr, N_TERM(ip_UsableSubsequences[k].q, i, j), true);
            }
         }
      }
   }
}

// Assign BJ64[i] = b^i (mod p) for each i in the ladder.
// Return b^Q (mod p).
uint32_t  CisOneWithMultipleSequencesWorker::ClimbLadder(MpArith mp, MpRes resBase)
{
   uint32_t  i, j, idx, lLen;

   lLen = *ip_AllLadders;
   
   // Precompute b^d (mod p) for 0 <= d <= Q, as necessary
   resBD[0] = mp.one();
   resBD[1] = resBase;
   resBD[2] = mp.mul(resBD[1], resBD[1]);
   
   i = 2;
   for (j=0; j<lLen; j++)
   {
      idx = ip_AllLadders[j+1];
      
      resBD[i+idx] = mp.mul(resBD[i], resBD[idx]);
      
      i += idx;
   }

   resBexpQ = resBD[ii_BestQ];

   return resBD[ii_BestQ];
}

// This function builds the list ii_CSSList[] of subsequences (k*b^d)*(b^Q)^m+c for
// which p may be a factor (-ckb^d is a quadratic/cubic/quartic/quintic
// residue with respect to p) and initialises the table D64[] with the
// values -c/(k*b^d) (mod p). As a side effect, bQ is set to the value b^Q
// (mod p) for use later in bsgs64(). Returns the number of subsequences listed in ii_CSS[].
uint32_t  CisOneWithMultipleSequencesWorker::SetupDiscreteLog(MpArith mp, MpRes resBase, MpRes resInvBase)
{
   uint64_t   bm, pShift, p = mp.p();
   uint32_t   idx, r;
   int32_t    shift;
   
   bm = p / 2;

   idx = bm % (ii_PowerResidueLcm/2);
   shift = ip_DivisorShifts[idx];
   
   if (shift == 0)
      return GetShift0Subsequences(mp, bm);
   
   if (shift > 0)
   {
      // p = 1 (mod s), where s is not a power of 2. Check for r-th power
      // residues for each prime power divisor r of s.
      pShift = p / shift;
   }
   else
   {
      // p = 1 (mod 2^s), where s > 1. Check for r-th power residues for each divisor
      // r of s. We handle this case seperately to avoid computing p/s using plain division.
      pShift = p >> (-shift);
      shift = 1 << (-shift);
   }

   /* For 0 <= r < s, resX[r] <- 1/(b^r)^((p-1)/s) */
   resX[0] = mp.one();
   resX[1] = mp.pow(resInvBase, pShift);
  
   for (r=1; resX[r] != resX[0]; r++)
      resX[r+1] = mp.mul(resX[r], resX[1]);

   if (shift % r != 0)
      FatalError("SetupDiscreteLog issue, shift %% xIdx != 0 (%u %% %u = %u)", shift, r, shift % r);
    
   // 1/(b^r)^((p-1)/s)=1 (mod p) therefore (1/(b^r)^((p-1)/s))^y=1 (mod p)
   // for 0 <= y < s/r. (Could we do more with this?)
   
   return GetShiftXSubsequences(mp, pShift, bm, r);
}

uint32_t  CisOneWithMultipleSequencesWorker::GetShift0Subsequences(MpArith mp, uint64_t bm)
{
   uint64_t   p = mp.p();
   int32_t    kcLegendre;
   bool       usable = false;
   uint32_t   qr_mod;
   uint32_t   j = 0;
   uint32_t   ssIdx;
   seq_t     *seqPtr;
   int32_t    bLegendre = 0;
   
   if (!ib_AllSequencesHaveLegendreTables)
      bLegendre = legendre(ii_Base, p);
      
   seqPtr = ip_FirstSequence;
   while (seqPtr != NULL)
   {
      legendre_t *legendrePtr = &ip_Legendre[seqPtr->seqIdx];
      
      if (legendrePtr->haveMap)
      {
         qr_mod = bm % legendrePtr->mod;
         usable = (legendrePtr->oneParityMap[L_BYTE(qr_mod)] & L_BIT(qr_mod));
      }
      else
      {
         kcLegendre = legendre(seqPtr->kcCore, p);
         
         switch (seqPtr->nParity)
         {
            case SP_EVEN:
               usable = (kcLegendre == 1);
               break;
            
            case SP_ODD:
               usable = (kcLegendre == bLegendre);
               break;
            
            case SP_MIXED:
               usable = (kcLegendre == 1 || kcLegendre == bLegendre);
               break;
            
            default:
               FatalError("parity not handled");
         }
      }
      
      if (usable)
      {
         uint64_t negCK = getNegCK(seqPtr, mp.p());
         MpRes resNegCK = mp.nToRes(negCK);
   
         for (ssIdx=seqPtr->ssIdxFirst; ssIdx<=seqPtr->ssIdxLast; ssIdx++)
         {
            ip_UsableSubsequences[j].seqPtr = seqPtr;
            ip_UsableSubsequences[j].q = ip_Subsequences[ssIdx].q;
            ip_UsableSubsequences[j].resBDCK = mp.mul(resBD[ip_Subsequences[ssIdx].q], resNegCK);

            j++;
         }
      }
      
      seqPtr = (seq_t *) seqPtr->next;
   }

   return j;
}

uint32_t  CisOneWithMultipleSequencesWorker::GetShiftXSubsequences(MpArith mp, uint64_t pShift, uint64_t bm, uint32_t r)
{
   uint64_t   p = mp.p();
   int32_t    kcLegendre;
   bool       usable = false;
   uint32_t   qr_mod, rIdx, idx;
   uint32_t   h, j = 0;
   uint32_t   ssIdx, cssIdx;
   seq_t     *seqPtr;
   int32_t    bLegendre = 0;
   
   rIdx = ip_PowerResidueIndices[r];
   
   if (!ib_AllSequencesHaveLegendreTables)
      bLegendre = legendre(ii_Base, p);

   seqPtr = ip_FirstSequence;
   while (seqPtr != NULL)
   {
      legendre_t *legendrePtr = &ip_Legendre[seqPtr->seqIdx];
      
      if (legendrePtr->haveMap)
      {
         qr_mod = bm % legendrePtr->mod;
         usable = (legendrePtr->oneParityMap[L_BYTE(qr_mod)] & L_BIT(qr_mod));
      }
      else
      {
         kcLegendre = legendre(seqPtr->kcCore, p);
         
         switch (seqPtr->nParity)
         {
            case SP_EVEN:
               usable = (kcLegendre == 1);
               break;
            
            case SP_ODD:
               usable = (kcLegendre == bLegendre);
               break;
            
            case SP_MIXED:
               usable = (kcLegendre == 1 || kcLegendre == bLegendre);
               break;
            
            default:
               FatalError("parity not handled");
         }
      }
      
      if (usable)
      {
         uint64_t negCK = getNegCK(seqPtr, mp.p());
         MpRes resNegCK = mp.nToRes(negCK);
         MpRes resPowNegCK = mp.pow(resNegCK, pShift);
           
         resX[r] = resPowNegCK;
         
         // Find h such that resX[h] = resX[r], i.e. (-ckb^h)^((p-1)/r)=1 (mod p), or h=r if not found
         for (h=0; resX[r] != resX[h]; h++)
            ;
    
         if (h < r)
         {
            qr_mod = cssIdx = CSS_INDEX(seqPtr->seqIdx, rIdx, h);

            // -c/(k*b^n) is an r-power residue for at least one term k*b^n+c of this sequence.
            cssIdx = ip_CongruentSubseqIndices[cssIdx];
  
  
            if (cssIdx > 0 && ip_AllSubseqs[cssIdx] > 0)
            {
               uint32_t *subseqs = &ip_AllSubseqs[cssIdx];
               uint32_t count = *subseqs;
               subseqs++;
     
               // -ckb^d is an r-th power residue for at least one term (k*b^d)*(b^Q)^(n/Q)+c of this subsequence.
               for (idx=0; idx<count; idx++)
               {
                  ssIdx = subseqs[idx];
                  
                  ip_UsableSubsequences[j].seqPtr = seqPtr;
                  ip_UsableSubsequences[j].q = ip_Subsequences[ssIdx].q;
                  ip_UsableSubsequences[j].resBDCK = mp.mul(resBD[ip_Subsequences[ssIdx].q], resNegCK);
                                    
                  j++;
               }
            }
         }
      }
      
      seqPtr = (seq_t *) seqPtr->next;
   }
   
   return j;
}

uint32_t  CisOneWithMultipleSequencesWorker::BabySteps(MpArith mp, MpRes resBase, MpRes resInvBase, uint32_t babySteps)
{
   MpRes    resInvBaseExpQ;
   MpRes    resBJ;
   MpRes    firstResBJ;
   uint32_t j;
      
   // b <- inv_b^Q (mod p)
   resInvBaseExpQ = mp.pow(resInvBase, ii_BestQ);
  
   firstResBJ = resBJ = mp.pow(resInvBaseExpQ, ii_SieveLow);

   for (j=0; j<babySteps; j++)
   { 
      ip_HashTable->Insert(resBJ, j);
      
      resBJ = mp.mul(resBJ, resInvBaseExpQ);
            
      if (resBJ == firstResBJ)
         return j + 1;
   }
   
   return 0;
}
