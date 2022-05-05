/* CisOneWithOneSequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

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
#include "CisOneWithOneSequenceHelper.h"
#include "CisOneWithOneSequenceWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "CisOneWithOneSequenceGpuWorker.h"
#endif

#define REPORT_STRFTIME_FORMAT "ETC %Y-%m-%d %H:%M"

CisOneWithOneSequenceHelper::CisOneWithOneSequenceHelper(App *theApp, uint64_t largestPrimeTested) : CisOneSequenceHelper(theApp, largestPrimeTested)
{
   theApp->WriteToConsole(COT_OTHER, "Sieving with single sequence c=1 logic for p >= %" PRIu64"", largestPrimeTested);

   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) theApp;
   
   if (srApp->GetBaseMultipleMulitplier() == 0)
      ii_BaseMultiple = 2 * DEFAULT_BM_MULTIPLIER_SINGLE;
   else
      ii_BaseMultiple = 2 * srApp->GetBaseMultipleMulitplier();
      
   if (srApp->GetPowerResidueLcmMultiplier() == 0)
      ii_PowerResidueLcm = ii_BaseMultiple * DEFAULT_PRL_MULTIPLIER_SINGLE;
   else
      ii_PowerResidueLcm = ii_BaseMultiple * srApp->GetPowerResidueLcmMultiplier();
   
   if (srApp->GetLimitBaseMultiplier() == 0)
      ii_LimitBase = ii_PowerResidueLcm * DEFAULT_LB_MULTIPLIER_SINGLE;
   else
      ii_LimitBase = ii_PowerResidueLcm * srApp->GetLimitBaseMultiplier();

   ip_App->WriteToConsole(COT_OTHER, "BASE_MULTIPLE = %u, POWER_RESIDUE_LCM = %u, LIMIT_BASE = %u", ii_BaseMultiple, ii_PowerResidueLcm, ii_LimitBase);
   
   ip_Legendre = NULL;
   ip_LegendreTable = NULL;

   BuildDivisorShifts();
   BuildPowerResidueIndices();
}

Worker  *CisOneWithOneSequenceHelper::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   AbstractWorker *theWorker;

   // Note that GenericWorker inherits from Worker.  This will not
   // only create the worker, but also start it.
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   if (gpuWorker)
      theWorker = new CisOneWithOneSequenceGpuWorker(id, ip_App, this);
   else
#endif
      theWorker = new CisOneWithOneSequenceWorker(id, ip_App, this);
      
   theWorker->Prepare(largestPrimeTested, ii_BestQ);

   return theWorker;
}

double   CisOneWithOneSequenceHelper::RateQ(uint32_t Q, uint32_t s)
{
   uint32_t          i;
   double            work;
   vector<double>    W;

   assert(Q % 2 == 0);
   assert(Q % ii_BaseMultiple == 0);
   assert(ii_LimitBase % Q == 0);

   W.resize(ii_PowerResidueLcm+1);

   work = 0;
   for (i = 2; i <= ii_PowerResidueLcm; i += 2)
   {
      if (ii_PowerResidueLcm % i == 0)
      {
         W[i] = EstimateWork(Q, (s+i-1)/i);

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

// These values originate from choose.c in sr1sieve.  The probably need to be
// adjusted for srsieve2, but I haven't put any thought into it.

// giantSteps are expensive compared to other loops, so it should
// have a height weight than the others
#define BABY_WORK    1.0    // 1 mulmod, 1 insert
#define GIANT_WORK   3.0    // 1 mulmod, 1 lookup
#define EXP_WORK     0.3    // 1 mulmod
#define SUBSEQ_WORK  1.4    // 1 mulmod, 1 lookup (giant step 0)
                               
// Q = q, s = number of subsequences.
double   CisOneWithOneSequenceHelper::EstimateWork(uint32_t Q, uint32_t s)
{
   uint32_t babySteps, giantSteps;
   double   work;

   ChooseSteps(Q, s, babySteps, giantSteps);

   work = babySteps*BABY_WORK + s*(giantSteps-1)*GIANT_WORK + Q*EXP_WORK + s*SUBSEQ_WORK;
   
   return work;
}

void  CisOneWithOneSequenceHelper::BuildCongruenceTables(void)
{
   uint64_t  bytesNeeded;
   seq_t    *seqPtr;
      
   ii_Dim3 = ii_PowerResidueLcm;
   ii_Dim2 = ii_Dim3 * ii_UsedPowerResidueIndices;
   ii_Dim1 = ii_Dim2 * SP_COUNT;

   ip_CongruentQIndices = (uint32_t *) xmalloc(ii_Dim1 * sizeof(uint32_t));
   ip_LadderIndices = (uint32_t *) xmalloc(ii_Dim1 * sizeof(uint32_t));

   // In sr1sieve, the congruent q and ladders are handled via four dimensional arrays.
   // For srsieve2, we will use two one dimensional arrays, which will be easier to pass to the GPU.
   // qIndices points to the first entry in qs for a specific parity, r, and h.
   // qs is a list of qs for that parity, r, and h with the first entry the length of the list
   // for that parity, r, and h.  The relationship for the ladder is the same.

   // The first entry if ip_AllQs and ip_AllLadders is unused.  This guarantees that a value of
   // 0 for the the qIndex or ladderIndex means that there are none. 
   ii_UsedQEntries = 1;
   ii_UsedLadderEntries = 1;
   ii_MaxQEntries = ii_MaxLadderEntries = ii_Dim1;
   
   // These will be resized as necessary
   ip_AllQs = (uint16_t *) xmalloc(ii_MaxQEntries * sizeof(uint16_t));
   ip_AllLadders = (uint16_t *) xmalloc(ii_MaxLadderEntries * sizeof(uint16_t));
   
   seqPtr = ip_FirstSequence;
   while (seqPtr != NULL)
   {
      BuildCongruenceTablesForSequence(seqPtr);
      
      seqPtr = (seq_t *) seqPtr->next;
   }
   
   bytesNeeded = ii_Dim1 * sizeof(uint32_t) * 2;
   ip_App->WriteToConsole(COT_OTHER, "%" PRIu64" bytes used for congruent q and ladder indices", bytesNeeded);
   
   bytesNeeded = (ii_MaxQEntries + ii_MaxLadderEntries) * sizeof(uint16_t);
   ip_App->WriteToConsole(COT_OTHER, "%" PRIu64" bytes used for congruent qs and ladders", bytesNeeded);
}

// Build tables sc_lists[i][r] of pointers to lists of subsequences
// whose terms k*b^m+c satisfy m = j (mod r)
void  CisOneWithOneSequenceHelper::BuildCongruenceTablesForSequence(seq_t *seqPtr)
{
   uint32_t   ssIdx, h, r, len[SP_COUNT];
   uint16_t  *tempQs[SP_COUNT];
   sp_t       parity = ip_FirstSequence->nParity;
   
   for (h=0; h<SP_COUNT; h++)
      tempQs[h] = (uint16_t *) xmalloc(ii_PowerResidueLcm * sizeof(uint16_t));
   
   for (r=1; r<=ii_PowerResidueLcm; r++)
   {
      if (ii_PowerResidueLcm % r != 0)
         continue;
      
      for (h=0; h<r; h++)
      {
         len[0] = len[1] = len[2] = 0;

         // TODO : can any subsequence be added twice?
         for (ssIdx=seqPtr->ssIdxFirst; ssIdx<=seqPtr->ssIdxLast; ssIdx++)
         {
            if (HasCongruentTerms(ssIdx, r, h))
            {
               uint16_t q = ip_Subsequences[ssIdx].q;
               
               // odd and even
               if (parity == SP_MIXED)
               {
                  tempQs[SP_MIXED][len[SP_MIXED]++] = q;
                     
                  if (ip_Subsequences[ssIdx].q%2 == 1)
                     tempQs[SP_ODD][len[SP_ODD]++] = q;

                  if (ip_Subsequences[ssIdx].q%2 == 0)
                     tempQs[SP_EVEN][len[SP_EVEN]++] = q;
               }
               
               if (parity == SP_ODD)
                  tempQs[SP_ODD][len[SP_ODD]++] = q;

               if (parity == SP_EVEN)
                  tempQs[SP_EVEN][len[SP_EVEN]++] = q;
            }
         }

         CopyQsAndMakeLadder(seqPtr, SP_EVEN,  r, h, tempQs[SP_EVEN],  len[SP_EVEN]);
         CopyQsAndMakeLadder(seqPtr, SP_ODD,   r, h, tempQs[SP_ODD],   len[SP_ODD]);
         CopyQsAndMakeLadder(seqPtr, SP_MIXED, r, h, tempQs[SP_MIXED], len[SP_MIXED]);
      }
   }
   
   for (h=0; h<SP_COUNT; h++)
      xfree(tempQs[h]);
}

void  CisOneWithOneSequenceHelper::CopyQsAndMakeLadder(seq_t *seqPtr, sp_t parity, uint32_t r, uint32_t h, uint16_t *qList, uint32_t qListLen)
{
   if (qListLen == 0)
      return;

   uint16_t rIdx = ip_PowerResidueIndices[r];
   uint32_t cqIdx = CQ_INDEX(parity, rIdx, h);
   
   ip_CongruentQIndices[cqIdx] = ii_UsedQEntries;
   ip_LadderIndices[cqIdx] = ii_UsedLadderEntries;

   // If the list of Qs isn't big enough, make it bigger
   if (ii_UsedQEntries + ii_BestQ + 10 >= ii_MaxQEntries)
   {
      uint32_t newMaxQEntries = ii_MaxQEntries + 1000;
      uint16_t *temp = (uint16_t *) xmalloc(newMaxQEntries * sizeof(uint16_t));

      memcpy(temp, ip_AllQs, ii_MaxQEntries * sizeof(uint16_t));
      xfree(ip_AllQs);

      ip_AllQs = temp;
      ii_MaxQEntries = newMaxQEntries;
   }

   // If the list of ladders isn't big enough, make it bigger
   if (ii_UsedLadderEntries + ii_BestQ + 10 >= ii_MaxLadderEntries)
   {
      uint32_t newMaxLadderEntries = ii_MaxLadderEntries + 1000;
      uint16_t *temp = (uint16_t *) xmalloc(newMaxLadderEntries * sizeof(uint16_t));

      memcpy(temp, ip_AllLadders, ii_MaxLadderEntries * sizeof(uint16_t));
      xfree(ip_AllLadders);

      ip_AllLadders = temp;
      ii_MaxLadderEntries = newMaxLadderEntries;
   }

   ip_AllQs[ii_UsedQEntries] = qListLen;
   ii_UsedQEntries++;
   
   memcpy(&ip_AllQs[ii_UsedQEntries], qList, qListLen * sizeof(uint16_t));
   ii_UsedQEntries += qListLen;
      
   MakeLadder(qList, qListLen);
}

// Note the the list of qs for the current sequence has already been filtered.
void   CisOneWithOneSequenceHelper::MakeLadder(uint16_t *qList, uint32_t qListLen)
{
   uint32_t          i, j, k, a;
   vector<uint8_t>   tempQs;
   
   assert(qListLen+1 < ii_BestQ);
   
   tempQs.resize(ii_BestQ+1);

   tempQs[ii_BestQ] = 1;
   
   for (i=0, a=1; i<qListLen; i++, a++)
      tempQs[qList[i]] = 1;

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

   a = 1;
   for (i=3; i<=ii_BestQ; i++)
      if (tempQs[i] == 2)
         a++;

   ip_AllLadders[ii_UsedLadderEntries] = a;
   ii_UsedLadderEntries++;
   
   j = 2;
   for (i=3; i<=ii_BestQ; i++)
      if (tempQs[i] == 2)
      {
         assert(tempQs[i-j]==2);
         ip_AllLadders[ii_UsedLadderEntries] = i - j;
         ii_UsedLadderEntries++;
         j = i;
      }
}

void  CisOneWithOneSequenceHelper::ComputeLegendreMemoryToAllocate(legendre_t *legendrePtr, uint64_t ssqfb)
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
            legendrePtr->bytesNeeded = (mapSize * 2);
         }
         break;
   }

   legendrePtr->r = r;
   legendrePtr->mod = mod;
   legendrePtr->canCreateMap = canCreateMap;
}

void     CisOneWithOneSequenceHelper::AssignMemoryToLegendreTable(legendre_t *legendrePtr, uint64_t bytesUsed)
{   
   switch (legendrePtr->nParity)
   {
      // odd n, test for (-bck/p)==1
      case SP_ODD: 
      case SP_EVEN:
         legendrePtr->oneParityMapIndex = bytesUsed;
         legendrePtr->oneParityMap = &ip_LegendreTable[bytesUsed];
         break;

      // odd and even n, test for (-ck/p)==1 and (-bck/p)==1
      default:
         legendrePtr->dualParityMapM1Index = bytesUsed;
         legendrePtr->dualParityMapM1 = &ip_LegendreTable[bytesUsed];

         bytesUsed += legendrePtr->mapSize;

         legendrePtr->dualParityMapP1Index = bytesUsed;
         legendrePtr->dualParityMapP1 = &ip_LegendreTable[bytesUsed];
         break;
   }
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
void  CisOneWithOneSequenceHelper::BuildLegendreTableForSequence(legendre_t *legendrePtr, uint64_t ssqfb, uint64_t stepsToDo, uint64_t stepsDone, time_t startTime)
{
   uint32_t i;
   uint32_t steps = legendrePtr->mod;
   int64_t  r = legendrePtr->r;
   double   percentDone;
   struct tm   *finish_tm;
   char     finishTimeBuffer[32];
   time_t   finishTime;
   
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
            if (jacobi(r, 2*i+1) == 1)
               legendrePtr->dualParityMapP1[L_BYTE(i)] |= L_BIT(i);
            
            if (jacobi(r*ssqfb, 2*i+1) == 1)
               legendrePtr->dualParityMapM1[L_BYTE(i)] |= L_BIT(i);
            
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
