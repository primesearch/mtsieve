/* CisOneWithMultipleSequencesHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include <time.h>
#include "../core/inline.h"
#include "SierpinskiRieselApp.h"
#include "CisOneWithMultipleSequencesHelper.h"
#include "CisOneWithMultipleSequencesWorker.h"

#ifdef HAVE_GPU_WORKERS
#include "CisOneWithOneSequenceGpuWorker.h"
#endif

#define REPORT_STRFTIME_FORMAT "ETC %Y-%m-%d %H:%M"

CisOneWithMultipleSequencesHelper::CisOneWithMultipleSequencesHelper(App *theApp, uint64_t largestPrimeTested) : CisOneSequenceHelper(theApp, largestPrimeTested)
{
   theApp->WriteToConsole(COT_OTHER, "Sieving with multi-sequence c=1 logic for p >= %" PRIu64"", largestPrimeTested);

   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) theApp;
   
   if (srApp->GetBaseMultipleMulitplier() == 0)
      ii_BaseMultiple = 2 * DEFAULT_BM_MULTIPLIER_MULTI;
   else
      ii_BaseMultiple = 2 * srApp->GetBaseMultipleMulitplier();
      
   if (srApp->GetPowerResidueLcmMultiplier() == 0)
      ii_PowerResidueLcm = ii_BaseMultiple * DEFAULT_PRL_MULTIPLIER_MULTI;
   else
      ii_PowerResidueLcm = ii_BaseMultiple * srApp->GetPowerResidueLcmMultiplier();
   
   if (srApp->GetLimitBaseMultiplier() == 0)
      ii_LimitBase = ii_PowerResidueLcm * DEFAULT_LB_MULTIPLIER_MULTI;
   else
      ii_LimitBase = ii_PowerResidueLcm * srApp->GetLimitBaseMultiplier();

   ip_App->WriteToConsole(COT_OTHER, "BASE_MULTIPLE = %u, POWER_RESIDUE_LCM = %u, LIMIT_BASE = %u", ii_BaseMultiple, ii_PowerResidueLcm, ii_LimitBase);
   
   ip_Legendre = NULL;
   ip_LegendreTable = NULL;

   BuildDivisorShifts();
   BuildPowerResidueIndices();
}

Worker  *CisOneWithMultipleSequencesHelper::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   AbstractWorker *theWorker;

   // Note that GenericWorker inherits from Worker.  This will not
   // only create the worker, but also start it.
   
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      theWorker = new CisOneWithOneSequenceGpuWorker(id, ip_App, this);
   else
#endif
      theWorker = new CisOneWithMultipleSequencesWorker(id, ip_App, this);
      
   theWorker->Prepare(largestPrimeTested, ii_BestQ);

   return theWorker;
}

double   CisOneWithMultipleSequencesHelper::RateQ(uint32_t Q, uint32_t s)
{
   vector<double>    W;
   double            work;
   uint32_t          i;
   //uint32_t r, babySteps, giantSteps;

   assert(Q % 2 == 0);
   assert(Q % ii_BaseMultiple == 0);
   assert(ii_LimitBase % Q == 0);

   W.resize(ii_PowerResidueLcm+1, 0);
   
   work = 0;
   for (i=2; i<=ii_PowerResidueLcm; i+=2)
   {
      if (ii_PowerResidueLcm % i == 0)
      {
         W[i] = EstimateWork(Q, s, i);

         // giantSteps are very expensive compared to other loops, so this might
         // always be the best choice even if it means more memory is needed
         //uint32_t r = 1 + ii_MaxN/Q - ii_MinN/Q;
         
         //ChooseSteps(r, s, babySteps, giantSteps);
         
         //W[i] = giantSteps;
      }

      if (gcd32(i+1, ii_PowerResidueLcm) == 1)
         work += W[gcd32(i, ii_PowerResidueLcm)];
   }

   return work;
}

// These values originate from choose.c in sr2sieve.  The probably need to be
// adjusted for srsieve2, but I haven't put any thought into it.

// giantSteps are expensive compared to other loops, so it should
// have a higher weight than the others
#define BABY_WORK    1.1    // 1 mulmod, 1 insert
#define GIANT_WORK   1.0    // 1 mulmod, 1 lookup
#define EXP_WORK     0.5    // 1 mulmod
#define SUBSEQ_WORK  1.0    // 1 mulmod, 1 lookup (giant step 0)
#define PRT_WORK     0.8    // List traversal, linear search, etc.
                               
// Q = q, s = number of subsequences.
double    CisOneWithMultipleSequencesHelper::EstimateWork(uint32_t Q, uint32_t s, uint32_t r)
{
   uint32_t babySteps, giantSteps;
   double   work;
   uint32_t x;

   x = (s+r-1)/r;

   ChooseSteps(Q, x, babySteps, giantSteps);

   work = babySteps*BABY_WORK + x*(giantSteps-1)*GIANT_WORK + Q*EXP_WORK + x*SUBSEQ_WORK;
   
   if (r > 2)
      work += x * PRT_WORK;
   
   return work;
}

void  CisOneWithMultipleSequencesHelper::BuildCongruenceTables(void)
{
   uint64_t  bytesNeeded;
   seq_t    *seqPtr;
   uint32_t  ssIdx, h, r, len;
   uint32_t *tempSubseqs;

   // Unlike the CisOneWithOneSequence, there is only one parity for this code
   ii_Dim3 = ii_PowerResidueLcm;
   ii_Dim2 = ii_Dim3 * ii_UsedPowerResidueIndices;
   ii_Dim1 = ii_Dim2 * ii_SequenceCount;

   ip_CongruentSubseqIndices = (uint32_t *) xmalloc(ii_Dim1 * sizeof(uint32_t));

   // In sr1sieve, the congruent q and ladders are handled via four dimensional arrays.
   // For srsieve2, we will use two one dimensional arrays, which will be easier to pass to the GPU.
   // qIndices points to the first entry in qs for a specific parity, r, and h.
   // qs is a list of qs for that parity, r, and h with the first entry the length of the list
   // for that parity, r, and h.  The relationship for the ladder is the same.
   
   ii_MaxSubseqEntries = 1000;
   ii_UsedSubseqEntries = 0;
   
   ip_AllSubseqs = (uint32_t *) xmalloc(ii_MaxSubseqEntries * sizeof(uint32_t));

   tempSubseqs = (uint32_t *) xmalloc(ii_PowerResidueLcm * sizeof(uint32_t));
      
   // Build tables sc_lists[i][r] of pointers to lists of subsequences
   // whose terms k*b^m+c satisfy m = j (mod r)
   for (r=1; r<=ii_PowerResidueLcm; r++)
   {
      if (ii_PowerResidueLcm % r != 0)
         continue;

      for (h=0; h<r; h++)
      {
         seqPtr = ip_FirstSequence;
         
         while (seqPtr != NULL)
         {
            len = 0;
            
            // TODO : can any subsequence be added twice?
            for (ssIdx=seqPtr->ssIdxFirst; ssIdx<=seqPtr->ssIdxLast; ssIdx++)
            {
               if (HasCongruentTerms(ssIdx, r, h))
               {
                  tempSubseqs[len] = ssIdx;
                  len++;
               }
            }

            CopySubseqs(seqPtr, r, h, tempSubseqs,  len);
            
            seqPtr = (seq_t *) seqPtr->next;
         }
      }
   }
   
   xfree(tempSubseqs);
   
   bytesNeeded = ii_Dim1 * sizeof(uint32_t);
   ip_App->WriteToConsole(COT_OTHER, "%" PRIu64" bytes used for congruent subseq indices", bytesNeeded);
   
   bytesNeeded = ii_MaxSubseqEntries * sizeof(uint32_t);
   ip_App->WriteToConsole(COT_OTHER, "%" PRIu64" bytes used for congruent subseqs", bytesNeeded);
   
   MakeLadder();
}

// In this class we are capturing the list of subsequences, not the list of qs.
// TODO:  AllSubseqs must be uin32_t for this, not uint16_t
void  CisOneWithMultipleSequencesHelper::CopySubseqs(seq_t *seqPtr, uint32_t r, uint32_t h, uint32_t *ssList, uint32_t ssListLen)
{
   if (ssListLen == 0)
      return;

   uint32_t rIdx = ip_PowerResidueIndices[r];
   uint32_t cssIdx = CSS_INDEX(seqPtr->seqIdx, rIdx, h);
      
   ip_CongruentSubseqIndices[cssIdx] = ii_UsedSubseqEntries;

   // If the list of Qs isn't big enough, make it bigger
   if (ii_UsedSubseqEntries + ssListLen + 10 >= ii_MaxSubseqEntries)
   {
      uint32_t newMaxSubseqEntries = ii_MaxSubseqEntries + 1000;
      uint32_t *temp = (uint32_t *) xmalloc(newMaxSubseqEntries * sizeof(uint32_t));

      memcpy(temp, ip_AllSubseqs, ii_MaxSubseqEntries * sizeof(uint32_t));
      xfree(ip_AllSubseqs);

      ip_AllSubseqs = temp;
      ii_MaxSubseqEntries = newMaxSubseqEntries;
   }
   
   ip_AllSubseqs[ii_UsedSubseqEntries] = ssListLen;
   ii_UsedSubseqEntries++;
   
   memcpy(&ip_AllSubseqs[ii_UsedSubseqEntries], ssList, ssListLen * sizeof(uint32_t));
   ii_UsedSubseqEntries += ssListLen;
}

void   CisOneWithMultipleSequencesHelper::MakeLadder()
{
   uint32_t          i, j, k, a, q;
   vector<uint8_t>   tempQs;
   
   tempQs.resize(ii_LimitBase, 0);
   tempQs[ii_BestQ] = 1;
   
   a = 1;
   for (i=0; i<ii_SubsequenceCount; i++)
   {     
      q = ip_Subsequences[i].q;
      if (tempQs[q] == 0)
      {
         tempQs[q] = 1;
         a++;
      }
   }
   
   for (i=0; i<3; i++)
   {
      if (tempQs[i] == 1)
        a--;
      tempQs[i] = 2;
   }

   while (a > 0)
   {
      for (i=3, j=2; i<=ii_BestQ; i++)
      {
         if (tempQs[i] == 2)
            j = i;
         else
            if (tempQs[i] == 1)
               break;
      }
      
      assert(i <= ii_BestQ);

      if (tempQs[i-j] == 2)
      {
         /* We can use an existing rung */
         tempQs[i] = 2;
         a--; 
      }
      else
      {
         /* Need to create a new rung */
         k = MIN(i-j,(i+1)/2); 
         assert(tempQs[k]==0);
         tempQs[k] = 1;
         a++;
         
         /* Need to re-check rungs above the new one */
         for (k++; k<=j; k++) 
            if (tempQs[k] == 2)
            {
               tempQs[k] = 1;
               a++;
            }
      }
   }

   a = 2;
   for (i=3; i<=ii_BestQ; i++)
      if (tempQs[i] == 2)
         a++;

   assert(a <= ii_BestQ);
   
   ip_AllLadders = (uint16_t *) xmalloc(a * sizeof(uint16_t));
   ii_LadderEntries = 1;
   
   j = 2;
   for (i=3; i<=ii_BestQ; i++)
      if (tempQs[i] == 2)
      {
         assert(tempQs[i-j]==2);
         ip_AllLadders[ii_LadderEntries] = i - j;
         ii_LadderEntries++;
         j = i;
      }
      
      
   ip_AllLadders[0] = ii_LadderEntries - 1;
}
 
void  CisOneWithMultipleSequencesHelper::ComputeLegendreMemoryToAllocate(legendre_t *legendrePtr, uint64_t ssqfb)
{
   int64_t     r = legendrePtr->kcCore;
   uint64_t    mod = legendrePtr->squareFreeK;
   uint32_t    mapSize = 0;
   bool        canCreateMap = true;
         
   switch (legendrePtr->nParity)
   {
      // odd n, test for (-bck/p)==1
      case SP_ODD: 
         mod *= ssqfb;
         r *= ssqfb;
         // Fall through

      // even n, test for (-ck/p)==1
      case SP_EVEN: 
         if ((r < 0 && (-r) % 4 != 3) || (r > 0 && r % 4 != 1))
            mod *= 2;
         
         mapSize = L_BYTES(mod);
               
         if (legendrePtr->squareFreeK > INT32_MAX/ssqfb)
            canCreateMap = false;
   
         if ((2*mod+1) >= INT32_MAX)
            canCreateMap = false;
         
         if (canCreateMap)
         {
            legendrePtr->mapSize = mapSize;
            legendrePtr->bytesNeeded = mapSize;
         }
         break;

      // odd and even n, test for (-ck/p)==1 and (-bck/p)==1
      default:
         mod = mod*2*ssqfb;
      
         mapSize = L_BYTES(mod);
         
         if (legendrePtr->squareFreeK > INT32_MAX/ssqfb)
            canCreateMap = false;
         
         if ((2*mod+1) >= INT32_MAX)
            canCreateMap = false;
         
         if (abs(r) >= INT32_MAX)
            canCreateMap = false;
         
         if (canCreateMap)
         {
            legendrePtr->mapSize = mapSize;
            legendrePtr->bytesNeeded = mapSize;
         }
         break;
   }
   
   legendrePtr->r = r;
   legendrePtr->mod = mod;
   legendrePtr->canCreateMap = canCreateMap;
}
 
void  CisOneWithMultipleSequencesHelper::AssignMemoryToLegendreTable(legendre_t *legendrePtr, uint64_t bytesUsed)
{
   legendrePtr->oneParityMapIndex = bytesUsed;
   legendrePtr->oneParityMap = &ip_LegendreTable[bytesUsed];
}

// For sequences with single parity terms (parity=+/-1):
// Set seq_mod and seq_map[0] for sequence k*b^n+c so that bit (p/2)%mod
// of seq_map[0] is set if and only if
// (-ck/p)=1 for sequences with all n even,
// (-bck/p)=1 for sequences with all n odd.
// 
// For sequences with mixed parity terms (parity=0):
// Set seq_mod, seq_map[0] and seq_map[1] for sequence k*b^n+c so that
// bit (p/2)%mod of seq_map[0] is set if and only if (-ck/p)=1,
// bit (p/2)%mod of seq_map[1] is set if and only if (-bck/p)=1.
// 
// In the worst case each table for k*b^n+c could be 4*b*k bits long.
void  CisOneWithMultipleSequencesHelper::BuildLegendreTableForSequence(legendre_t *legendrePtr, uint64_t ssqfb, uint64_t stepsToDo, uint64_t stepsDone, time_t startTime)
{
   uint32_t    i;
   uint32_t    steps = legendrePtr->mod;
   int64_t     r = legendrePtr->r;
   double      percentDone;
   struct tm  *finish_tm;
   char        finishTimeBuffer[32];
   time_t      finishTime;
      
   switch (legendrePtr->nParity)
   {
      // odd n, test for (-bck/p)==1
      case SP_ODD: 
      case SP_EVEN:        
         for (i=0; i<steps; i++)
         {
            if (jacobi(r, 2*i+1) == 1)
               legendrePtr->oneParityMap[L_BYTE(i)] |= L_BIT(i);

            if (!(++stepsDone & 0xffffff))
            {
               percentDone = ((double) stepsDone)/stepsToDo;
               finishTime = (time_t) (startTime + (time(NULL)-startTime)/percentDone);
               
               finish_tm = localtime(&finishTime);
               if (!finish_tm || !strftime(finishTimeBuffer, sizeof(finishTimeBuffer), REPORT_STRFTIME_FORMAT, finish_tm))
                  finishTimeBuffer[0] = '\0';
      
               ip_App->WriteToConsole(COT_SIEVE, "Building Legendre tables: %.1f%% done %s (currently at k=%" PRIu64")", 100.0*percentDone, finishTimeBuffer, legendrePtr->k);
            }
         }
         break;

      // odd and even n, test for (-ck/p)==1 and (-bck/p)==1
      default:         
         for (i=0; i<steps; i++)
         {
            if (jacobi(r, 2*i+1) == 1 || jacobi(r*ssqfb, 2*i+1) == 1)
               legendrePtr->oneParityMap[L_BYTE(i)] |= L_BIT(i);

            if (!(++stepsDone & 0xffffff))
            {
               percentDone = ((double) stepsDone)/stepsToDo;
               finishTime = (time_t) (startTime + (time(NULL)-startTime)/percentDone);
               
               finish_tm = localtime(&finishTime);
               if (!finish_tm || !strftime(finishTimeBuffer, sizeof(finishTimeBuffer), REPORT_STRFTIME_FORMAT, finish_tm))
                  finishTimeBuffer[0] = '\0';
      
               ip_App->WriteToConsole(COT_SIEVE, "Building Legendre tables: %.1f%% done %s (currently at k=%" PRIu64")", 100.0*percentDone, finishTimeBuffer, legendrePtr->k);
            }
         }
         break;
   }
}
