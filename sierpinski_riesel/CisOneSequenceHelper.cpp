/* CisOneSequenceHelper.cpp -- (C) Mark Rodenkirch, January 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include "CisOneSequenceHelper.h"
#include "CisOneWithOneSequenceWorker.h"
#include "../core/inline.h"

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

CisOneSequenceHelper::CisOneSequenceHelper(App *theApp, uint64_t largestPrimeTested) : AbstractSequenceHelper(theApp, largestPrimeTested)
{
   // The parent's constuctor does much of the work

   uint32_t shifts = POWER_RESIDUE_LCM / 2;
   uint32_t i, r, idx;
   int32_t  shift, divide;
   
   ip_DivisorShifts = (int16_t *) xmalloc(shifts * sizeof(int16_t));

   for (i=0; i<shifts; i++)
   {
      r = gcd32(2*i, POWER_RESIDUE_LCM);
      
      divide = r;
      for (shift = 0; divide % 2 == 0; shift++)
         divide /= 2;
         
      if (divide > 1)
      {
         ip_DivisorShifts[i] = r;
         continue;
      }

      if (shift > 1)
         ip_DivisorShifts[i] = -shift;
      else
         ip_DivisorShifts[i] = 0;
   }
   
   ii_PrlCount = 0;
   ii_LadderCount = 0;
   ii_MaxLadderRungs = 0;
   ii_MaxQs = 0;

   for (r=1; r<=POWER_RESIDUE_LCM; r++)
   {
      if (POWER_RESIDUE_LCM % r != 0)
         continue;
         
      ii_PrlCount++;
   }
   
   ip_PrlIndices = (uint16_t *) xmalloc((POWER_RESIDUE_LCM+1) * sizeof(uint16_t));
   idx = 1;

   for (r=1; r<=POWER_RESIDUE_LCM; r++)
   {
      if (POWER_RESIDUE_LCM % r != 0)
         continue;

      ip_PrlIndices[r] = idx;
      idx++;
   }         
}

void   CisOneSequenceHelper::CleanUp(void)
{
   seq_t    *seq;

   xfree(ip_PrlIndices);
   
   seq = ip_FirstSequence;
   while (seq != NULL)
   {     
      if (seq->congruentQs != NULL)
         xfree(seq->congruentQs);
      
      if (seq->congruentQIndices != NULL)
         xfree(seq->congruentQIndices);
         
      if (seq->congruentLadders != NULL)
         xfree(seq->congruentLadders);
 
      if (seq->congruentLadderIndices != NULL)
         xfree(seq->congruentLadderIndices);
 
      seq = (seq_t *) seq->next;
   }
}

void        CisOneSequenceHelper::LastChanceLogicBeforeSieving(void)
{
   uint32_t  ssIdx, babySteps, giantSteps;
   uint32_t  sieveLow = ii_MinN / ii_BestQ;
   uint32_t  sieveHigh = ii_MaxN / ii_BestQ;
   uint64_t  preBuildCpuBytes, postBuildCpuBytes;
   seq_t    *seq;
         
   ii_MaxBabySteps = 0;
   
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      ChooseSteps(ii_BestQ, ssIdx+1, babySteps, giantSteps);
      
      if (sieveHigh >= sieveLow + babySteps*giantSteps)
         FatalError("LastChanceLogicBeforeSieving miscomputed the steps");
         
      ip_Subsequences[ssIdx].babySteps = babySteps;
      ip_Subsequences[ssIdx].giantSteps = giantSteps;

      if (ii_MaxBabySteps < babySteps)
         ii_MaxBabySteps = babySteps;
   }
   
   // These are much bigger than we need them to be.  We will shrink and copy to
   // seq->congruentQs and seq->congruentLadders when done
   ip_TempQs = (uint16_t *) xmalloc(100000000 * sizeof(uint16_t));
   ip_TempLadders = (uint16_t *) xmalloc(100000000 * sizeof(uint16_t));
   
   preBuildCpuBytes = GetCpuMemoryUsage();
   
   seq = ip_FirstSequence;
   while (seq != NULL)
   {     
      MakeSubseqCongruenceTables(seq);
      
      seq = (seq_t *) seq->next;
   }
   
   postBuildCpuBytes = GetCpuMemoryUsage();

   xfree(ip_TempQs);
   xfree(ip_TempLadders);  
   
   ip_App->WriteToConsole(COT_OTHER, "%llu bytes used for congruence tables", postBuildCpuBytes - preBuildCpuBytes);
      
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   
   ib_UseLegendreTables = srApp->UseLegendreTables();

   if (ib_UseLegendreTables)
      ib_UseLegendreTables = BuildLegendreTables(srApp->GetLegendreFileName());
}

uint32_t    CisOneSequenceHelper::FindBestQ(uint32_t &expectedSubsequences)
{
   uint32_t   i = 0, j, n;
   uint32_t   bit;
   uint32_t   S[NDIVISORS];
   uint32_t   W[NDIVISORS];
   bool       R[LIMIT_BASE];
   seq_t     *seq;

   seq = ip_FirstSequence;
   do
   {
      for (j=0; j<NDIVISORS; j++)
         S[j] = 0;

      for (j=0; j<LIMIT_BASE; j++)
         R[j] = false;
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[bit])
            R[n%LIMIT_BASE] = true;
         
         bit++;
      }
 
      i = 0;
      for (j=0; j<NDIVISORS; j++)
      {
         if (NDIVISORS % (j+1) == 0)
         {
            S[j] += CountResidueClasses((j+1)*BASE_MULTIPLE, LIMIT_BASE, R);
            
            W[j] = RateQ((j+1)*BASE_MULTIPLE, S[j]);
            
            if (W[j] < W[i])
               i = j;
         }
      }

      seq = (seq_t *) seq->next;
   } while (seq != NULL);

   expectedSubsequences = S[i];
   
   return (i+1)*BASE_MULTIPLE;
}

bool   CisOneSequenceHelper::BuildLegendreTables(string legendreFileName)
{
   uint64_t     preBuildCpuBytes;
   uint64_t     postBuildCpuBytes;
   legendre_t  *legendre;
   bool         builtLegendreTables;
   seq_t       *seq;
      
   preBuildCpuBytes = GetCpuMemoryUsage();
   
   seq = ip_FirstSequence;
   do
   {
      legendre = (legendre_t *) xmalloc(sizeof(legendre_t));
      
      seq->legendrePtr = legendre;

      builtLegendreTables = BuildLegendreTableForSequence(seq, legendre);
      
      if (!builtLegendreTables)
         break;

      seq = (seq_t *) seq->next;
   } while (seq != NULL);

   postBuildCpuBytes = GetCpuMemoryUsage();
   
   if (builtLegendreTables)
      ip_App->WriteToConsole(COT_OTHER, "%llu bytes used for Legendre tables", postBuildCpuBytes - preBuildCpuBytes);

   return builtLegendreTables;
}
