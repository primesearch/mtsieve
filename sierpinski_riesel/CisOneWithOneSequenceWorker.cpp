/* CisOneWithOneSequenceWorker.cpp -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "CisOneWithOneSequenceWorker.h"
#include "../core/inline.h"
#include "../core/MpArith.h"

#define N_TERM(q, i, j)      ((ii_SieveLow + (j) + (i)*babySteps)*ii_BestQ + q)

CisOneWithOneSequenceWorker::CisOneWithOneSequenceWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ip_CisOneHelper = (CisOneWithOneSequenceHelper *) appHelper;
   
   // Everything we need is done in the constuctor of the parent class
   ib_Initialized = true;

   ii_BaseMultiple = ip_CisOneHelper->GetBaseMultiple();
   ii_LimitBase = ip_CisOneHelper->GetLimitBase();
   ii_PowerResidueLcm = ip_CisOneHelper->GetPowerResidueLcm();
}

void  CisOneWithOneSequenceWorker::CleanUp(void)
{
   delete ip_HashTable;
   
   xfree(resBD);
   xfree(resBJ);
   xfree(resX);
}

void  CisOneWithOneSequenceWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{
   ii_BestQ = bestQ;   
   ii_SieveLow = ii_MinN / ii_BestQ;

   resX = (MpRes *) xmalloc((ii_PowerResidueLcm+5) * sizeof(MpRes));
   resBD = (MpRes *) xmalloc((ii_SubsequenceCount+4)*sizeof(MpRes));
   resBJ = (MpRes *) xmalloc((ii_SubsequenceCount+4)*sizeof(MpRes));
   
   ip_DivisorShifts = ip_CisOneHelper->GetDivisorShifts();
   ip_PowerResidueIndices = ip_CisOneHelper->GetPowerResidueIndices();
      
   ii_Dim1 = ip_CisOneHelper->GetDim1();
   ii_Dim2 = ip_CisOneHelper->GetDim2();
   ii_Dim3 = ip_CisOneHelper->GetDim3();
   
   ip_CongruentQIndices = ip_CisOneHelper->GetCongruentQIndices();
   ip_LadderIndices = ip_CisOneHelper->GetLadderIndices();
   
   ip_AllQs = ip_CisOneHelper->GetAllQs();
   ip_AllLadders = ip_CisOneHelper->GetAllLadders();
   
   ip_Legendre = ip_CisOneHelper->GetLegendre();
   ip_LegendreTable = ip_CisOneHelper->GetLegendreTable();
   
   ip_HashTable = new HashTable(ip_CisOneHelper->GetMaxBabySteps());
}

void  CisOneWithOneSequenceWorker::TestMegaPrimeChunk(void)
{
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint64_t p;
   sp_t     parity;

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

      parity = GetParity(p);
      
      if (parity != SP_NO_PARITY)
         TestSinglePrime(p, parity);

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
void  CisOneWithOneSequenceWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CisOneWithOneSequenceWorker::TestMiniPrimeChunk not implemented");
}

sp_t   CisOneWithOneSequenceWorker::GetParity(uint64_t p)
{
   legendre_t *legendrePtr = ip_Legendre;
      
   // Mixed parity sequences
   if (ip_FirstSequence->nParity == SP_MIXED)
   {
      bool qr_m1, qr_p1;
      
      if (legendrePtr->haveMap)
      {
         uint32_t qr_mod = (p/2) % legendrePtr->mod;
                  
         qr_m1 = (legendrePtr->dualParityMapM1[L_BYTE(qr_mod)] & L_BIT(qr_mod));
         qr_p1 = (legendrePtr->dualParityMapP1[L_BYTE(qr_mod)] & L_BIT(qr_mod));
      }
      else
      {
         int32_t sym = legendre(ip_FirstSequence->kcCore, p);
         
         qr_p1 = (sym == 1);
         qr_m1 = (sym == legendre(ii_Base, p));
      }
      
      if (qr_m1)
         return (qr_p1 ? SP_MIXED : SP_ODD);

      return (qr_p1 ? SP_EVEN : SP_NO_PARITY);
   }
   
   // Single parity sequences
   bool qr;
   
   if (legendrePtr->haveMap)
   {
      uint32_t qr_mod = (p/2) % legendrePtr->mod;
         
      qr = (legendrePtr->oneParityMap[L_BYTE(qr_mod)] & L_BIT(qr_mod));
   }
   else
   {
      int32_t sym = legendre(ip_FirstSequence->kcCore, p);
      
      if (ip_FirstSequence->nParity == SP_EVEN)
         qr = (sym == 1);
      else
         qr = (sym == legendre(ii_Base, p));
   }
   
   if (qr)
      return ip_FirstSequence->nParity;
      
   return SP_NO_PARITY;
}

void  CisOneWithOneSequenceWorker::TestSinglePrime(uint64_t p, sp_t parity)
{
   uint64_t   invBase, negCK;
   uint32_t   k, orderOfB, ssCount;
   uint32_t   babySteps, giantSteps;
   uint32_t   i, j;
   uint32_t   cqIdx, qIdx;
   uint16_t  *seqQs;
   
   MpArith mp(p);

   // compute 1/base (mod p)
   invBase = invmod64(ii_Base, p);

   negCK = getNegCK(ip_FirstSequence, p);
      
   MpRes resBase = mp.nToRes(ii_Base);
   MpRes resInvBase = mp.nToRes(invBase);
   MpRes resNegCK = mp.nToRes(negCK);
      
   cqIdx = SetupDiscreteLog(mp, resBase, resInvBase, resNegCK, parity);
   
   qIdx = ip_CongruentQIndices[cqIdx];
      
   // If no qs for this p, then no factors, so return
   if (qIdx == 0)
      return;

   seqQs = &ip_AllQs[qIdx];
      
   // The number of subsequences (qs) for this sequence
   ssCount = seqQs[0];
   
   // If no subsequences for this p, then no factors, so return
   if (ssCount == 0)
      return;
   
   // Skip the count since we have already copied it
   seqQs++;
      
   // -ckb^d is an r-th power residue for at least one term (k*b^d)*(b^Q)^(n/Q)+c of this subsequence
   BuildLookupsAndClimbLadder(mp, resBase, resNegCK, cqIdx, ssCount, seqQs);
   
   ip_HashTable->Clear();

   babySteps = ip_Subsequences[ssCount-1].babySteps;
   giantSteps = ip_Subsequences[ssCount-1].giantSteps;

   orderOfB = BabySteps(mp, resBase, resInvBase, babySteps);
   
   if (orderOfB > 0)
   {
      // If orderOfB > 0, then this is all the information we need to
      // determine every solution for this p, so no giant steps are neede
      for (k=0; k<ssCount; k++)
      {
          j = ip_HashTable->Lookup(resBD[k]);

          while (j < babySteps * giantSteps)
          {
             ip_SierpinskiRieselApp->ReportFactor(p, ip_FirstSequence, N_TERM(seqQs[k], 0, j), true);
             
             j += orderOfB;
          }
      }
   }
   else
   {
      // First giant step
      for (k=0; k<ssCount; k++)
      {
         j = ip_HashTable->Lookup(resBD[k]);

         if (j != HASH_NOT_FOUND)
            ip_SierpinskiRieselApp->ReportFactor(p, ip_FirstSequence, N_TERM(seqQs[k], 0, j), true);
      }

      // Remaining giant steps
      if (giantSteps > 1)
      {
         MpRes resBQM = mp.pow(resBexpQ, babySteps);
         
         for (i=1; i<giantSteps; i++)
         {
            for (k=0; k<ssCount; k++)
            {
               resBD[k] = mp.mul(resBD[k], resBQM);
               
               j = ip_HashTable->Lookup(resBD[k]);

               if (j != HASH_NOT_FOUND)
                  ip_SierpinskiRieselApp->ReportFactor(p, ip_FirstSequence, N_TERM(seqQs[k], i, j), true);
            }
         }
      }
   }
}

// This function builds the list ii_CSSList[] of subsequences (k*b^d)*(b^Q)^m+c for
// which p may be a factor (-ckb^d is a quadratic/cubic/quartic/quintic
// residue with respect to p) and initialises the table D64[] with the
// values -c/(k*b^d) (mod p). As a side effect, bQ is set to the value b^Q
// (mod p) for use later in bsgs64(). Returns the number of subsequences listed in ii_CSS[].
uint32_t  CisOneWithOneSequenceWorker::SetupDiscreteLog(MpArith mp, MpRes resBase, MpRes resInvBase, MpRes resNegCK, sp_t parity)
{
   uint64_t   pShift, p = mp.p();
   uint32_t   idx;
   uint32_t   h, r;
   uint32_t   rIdx;
   int16_t    shift;
   
   idx = (p/2) % (ii_PowerResidueLcm/2);
   shift = ip_DivisorShifts[idx];

   if (shift == 0)
   {
      rIdx = ip_PowerResidueIndices[1];
   
      return CQ_INDEX(parity, rIdx, 0);
   }
   
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
     
   resX[0] = mp.one();
   resX[1] = mp.pow(resInvBase, pShift);
  
   for (r=1; resX[r] != resX[0]; r++)
      resX[r+1] = mp.mul(resX[r], resX[1]);

   if (shift % r != 0)
      FatalError("SetupDiscreteLog issue, shift %% xIdx != 0 (%u %% %u = %u)", shift, r, shift % r);
    
   // resX[r] <- (-ck)^((p-1)/shift)
   resX[r] = mp.pow(resNegCK, pShift);

   // Find h such that resX[h] = resX[r], i.e. (-ckb^h)^((p-1)/r)=1 (mod p), or h=r if not found
   for (h=0; resX[r] != resX[h]; h++)
      ;

   // If no h was found, then there is nothing further to do.
   if (h == r)
      return 0;

   rIdx = ip_PowerResidueIndices[r];
   
   return CQ_INDEX(parity, rIdx, h);
}

// Assign BJ64[i] = b^i (mod p) for each i in the ladder.
// Return b^Q (mod p).
void  CisOneWithOneSequenceWorker::BuildLookupsAndClimbLadder(MpArith mp, MpRes resBase, MpRes resNegCK, uint16_t cqIdx, uint16_t ssCount, uint16_t *seqQs)
{
   uint32_t  i, j, idx, lLen, ladderIdx;
   uint16_t *ladders;

   ladderIdx = ip_LadderIndices[cqIdx];
   ladders = &ip_AllLadders[ladderIdx];
   
   lLen = ladders[0];
   ladders++;
   
   // Precompute b^d (mod p) for 0 <= d <= Q, as necessary
   resX[0] = mp.one();
   resX[1] = resBase;
   resX[2] = mp.mul(resX[1], resX[1]);
   
   i = 2;
   for (j=0; j<lLen; j++)
   {
      idx = ladders[j];
      
      resX[i+idx] = mp.mul(resX[i], resX[idx]);
      
      i += idx;
   }
   
   resBexpQ = resX[ii_BestQ];

   for (j=0; j<ssCount; j++)
      resBD[j] = mp.mul(resX[seqQs[j]], resNegCK);
}

uint32_t  CisOneWithOneSequenceWorker::BabySteps(MpArith mp, MpRes resBase, MpRes resInvBase, uint32_t babySteps)
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
