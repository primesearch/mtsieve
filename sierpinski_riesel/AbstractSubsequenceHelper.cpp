/* AbstractSubsequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "SierpinskiRieselApp.h"
#include "AbstractSubsequenceHelper.h"

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

AbstractSubsequenceHelper::AbstractSubsequenceHelper(App *theApp, seq_t *sequences, uint32_t sequenceCount)
{
   ip_App = theApp;
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) theApp;
   
   ii_Base = srApp->GetBase();
   ii_MinN = srApp->GetMinN();
   ii_MaxN = srApp->GetMaxN();
   
   ii_BestQ = 1;
   ii_MinM = srApp->GetMinN();
   ii_MaxM = srApp->GetMaxN();

   ip_Sequences = sequences;
   ii_SequenceCount = sequenceCount;
   
   ip_Subsequences = 0;
   ii_SubsequenceCount = 0;
   ii_SubsequenceCapacity = 0;
}

void     AbstractSubsequenceHelper::CleanUp(void)
{
   xfree(ip_Subsequences);
}

uint64_t    AbstractSubsequenceHelper::MakeSubsequencesForNewSieve(seq_t *sequences)
{
   uint32_t seqIdx;
   uint32_t ssIdx;
   uint32_t mTerms = (ii_MaxN - ii_MinN + 1);
   uint64_t termCount = 0;

   // For a new sieve there is only one subsequence per sequence
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      ssIdx = AddSubsequence(seqIdx, 0, mTerms);
      
      subseq_t *ssPtr = &ip_Subsequences[ssIdx];
 
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (ip_Sequences[seqIdx].nTerms[NBIT(n)])
         {
            ssPtr->mTerms[NBIT(n)] = true;
            termCount++;
         }
      }
   }

   ii_BestQ = 1;
   
   return termCount;
}

void      AbstractSubsequenceHelper::CreateEmptySubsequences(uint32_t subsequenceCount)
{
   ip_Subsequences = (subseq_t *) xmalloc(subsequenceCount * sizeof(subseq_t));
   
   ii_SubsequenceCount = 0;
   ii_SubsequenceCapacity = subsequenceCount;
}

uint32_t  AbstractSubsequenceHelper::AddSubsequence(uint32_t seqIdx, uint32_t d, uint32_t mTermCount)
{
   uint32_t ssIdx;
   
   ip_Sequences[seqIdx].ssCount++;
   
   ssIdx = ii_SubsequenceCount;
   ii_SubsequenceCount++;
   
   if (ip_Sequences[seqIdx].ssCount == 1)
   {
      ip_Sequences[seqIdx].ssIdxFirst = ii_SubsequenceCount;
      ip_Sequences[seqIdx].parity = d % 2;
   }
   
   if (ip_Sequences[seqIdx].parity != d % 2)
      ip_Sequences[seqIdx].parity = 2;
   
   ip_Sequences[seqIdx].ssIdxLast = ii_SubsequenceCount;
   
   subseq_t *ssPtr = &ip_Subsequences[ssIdx];
   
   ssPtr->seqIdx = seqIdx;
   ssPtr->k = ip_Sequences[seqIdx].k;
   ssPtr->c = ip_Sequences[seqIdx].c;
   ssPtr->d = d;
   ssPtr->mTerms.resize(mTermCount);
      
   std::fill(ssPtr->mTerms.begin(), ssPtr->mTerms.end(), false);

   return ssIdx;
}
