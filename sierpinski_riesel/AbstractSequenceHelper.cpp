/* AbstractSequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include <cfenv>
#include "SierpinskiRieselApp.h"
#include "AbstractSequenceHelper.h"
#include "../core/inline.h"
#include "../core/SmallHashTable.h"

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

AbstractSequenceHelper::AbstractSequenceHelper(App *theApp, uint64_t largestPrimeTested)
{
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) theApp;
   
   ip_App = (SierpinskiRieselApp *) theApp;
   
   ii_Base = srApp->GetBase();
   ii_MinN = srApp->GetMinN();
   ii_MaxN = srApp->GetMaxN();
   
   ii_BestQ = 1;
   ii_MinM = srApp->GetMinN();
   ii_MaxM = srApp->GetMaxN();
   
   ip_FirstSequence = srApp->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   
   ip_Subsequences = 0;
   ii_SubsequenceCount = 0;
   ii_SubsequenceCapacity = 0;
   ii_MaxSubsequenceCount = 0;
}

void     AbstractSequenceHelper::CleanUp(void)
{
   xfree(ip_Subsequences);
}

uint64_t    AbstractSequenceHelper::MakeSubsequencesForNewSieve(void)
{
   uint32_t    ssIdx;
   uint32_t    nTerms = (ii_MaxN - ii_MinN + 1);
   uint64_t    termCount = 0;
   seq_t      *seqPtr;

   CreateEmptySubsequences(ii_SequenceCount);
   
   seqPtr = ip_FirstSequence;
   do
   {
      ssIdx = AddSubsequence(seqPtr, 0, nTerms);
      
      subseq_t *ssPtr = &ip_Subsequences[ssIdx];
 
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seqPtr->nTerms[NBIT(n)])
         {
            ssPtr->mTerms[NBIT(n)] = true;
            termCount++;
         }
      }
      
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);

   ii_BestQ = 1;
   
   return termCount;
}

void        AbstractSequenceHelper::MakeSubsequencesForOldSieve(uint64_t expectedTerms)
{
   std::vector<bool>     needss; 
   std::vector<uint32_t> rss;
   
   uint32_t  bit, r, n;
   uint32_t  expectedSubsequences;
   uint64_t  countedTerms = 0;
   uint32_t  ssIdx;
   seq_t    *seqPtr;
   
   if (ip_Subsequences)
      xfree(ip_Subsequences);
   
   ip_Subsequences = 0;
   ii_SubsequenceCount = 0;
   ii_SubsequenceCapacity = 0;
   ii_MaxSubsequenceCount = 0;
   
   ii_BestQ = FindBestQ(expectedSubsequences);
   
   ii_MinM = ii_MinN/ii_BestQ;
   ii_MaxM = ii_MaxN/ii_BestQ;

   needss.resize(ii_BestQ);
   rss.resize(ii_BestQ);

   // This is the maximum number of subsequences
   CreateEmptySubsequences(ii_SequenceCount * ii_BestQ);
   
   seqPtr = ip_FirstSequence;
   do
   {
      seqPtr->ssCount = 0;
      seqPtr->ssIdxFirst = 0;
      seqPtr->ssIdxLast = 0;
      
      std::fill(needss.begin(), needss.end(), false);

      bit = NBIT(ii_MinN);
      
      // Determine which subsequences we need to build
      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seqPtr->nTerms[bit])
         {
            r = n % ii_BestQ;
            needss[r] = true;
         }
         
         bit++;
      }

      for (r=0; r<ii_BestQ; r++)
      {
         if (needss[r])
         {
            rss[r] = ii_SubsequenceCount;
            AddSubsequence(seqPtr, r, ii_MaxM - ii_MinM + 1);
         }
      }
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seqPtr->nTerms[bit])
         {
            r = n % ii_BestQ;
            ssIdx = rss[r];
            ip_Subsequences[ssIdx].mTerms[n/ii_BestQ - ii_MinM] = true;
            countedTerms++;
         }
         
         bit++;
      }
      
      seqPtr = (seq_t *) seqPtr->next;
   } while (seqPtr != NULL);
   
   if (ii_SubsequenceCount != expectedSubsequences)
      FatalError("Expected %u subsequences but %u were created", expectedSubsequences, ii_SubsequenceCount);

   if (expectedTerms != countedTerms)
      FatalError("Expected %" PRIu64" terms when building sequences, but counted only %" PRIu64"", expectedTerms, countedTerms);
   
   const char *sequenceText = ((ii_SequenceCount > 1) ? "sequences" : "sequence");
  
   ip_App->WriteToConsole(COT_OTHER, "Split %u base %u %s into %u base %u^%u sequences.", ii_SequenceCount,
		ii_Base, sequenceText, ii_SubsequenceCount, ii_Base, ii_BestQ);
}

void      AbstractSequenceHelper::CreateEmptySubsequences(uint32_t subsequenceCount)
{
   ip_Subsequences = (subseq_t *) xmalloc(subsequenceCount, sizeof(subseq_t), "subsequences");
   
   ii_SubsequenceCount = 0;
   ii_SubsequenceCapacity = subsequenceCount;
}

uint32_t  AbstractSequenceHelper::AddSubsequence(seq_t *seqPtr, uint32_t q, uint32_t mTermCount)
{
   uint32_t ssIdx;

   seqPtr->ssCount++;
   
   if (seqPtr->ssCount == 1)
   {
      seqPtr->ssIdxFirst = ii_SubsequenceCount;
      
      if (q % 2 == 0)
         seqPtr->nParity = SP_EVEN;
      else
         seqPtr->nParity = SP_ODD;
   }
   else{
      if (seqPtr->nParity == SP_EVEN && q % 2 != 0)
         seqPtr->nParity = SP_MIXED;
      
      if (seqPtr->nParity == SP_ODD && q % 2 == 0)
         seqPtr->nParity = SP_MIXED;
   }
   
   // This is the index of the new subsequence being added
   ssIdx = ii_SubsequenceCount;

   seqPtr->ssIdxLast = ssIdx;
   
   subseq_t *ssPtr = &ip_Subsequences[ssIdx];
   
   ssPtr->seqPtr = seqPtr;
   ssPtr->k = seqPtr->k;
   ssPtr->c = seqPtr->c;
   ssPtr->q = q;
   ssPtr->mTerms.resize(mTermCount, false);

   if (seqPtr->ssCount > ii_MaxSubsequenceCount)
      ii_MaxSubsequenceCount = seqPtr->ssCount;

   ii_SubsequenceCount++;
   return ssIdx;
}

uint32_t    AbstractSequenceHelper::CountResidueClasses(uint32_t d, uint32_t Q, std::vector<bool> R)
{
   uint32_t i, count;
   std::vector<bool>  R0;

   assert(Q % d == 0);

   R0.resize(d, false);

   for (i = 0; i < Q; i++)
      if (R[i])
         R0[i%d] = true;

   for (i = 0, count = 0; i < d; i++)
      if (R0[i])
         count++;
   
   return count;
}

void  AbstractSequenceHelper::ChooseSteps(uint32_t Q, uint32_t s, uint32_t &babySteps, uint32_t &giantSteps)
{
   uint32_t m, M;
   int32_t  roundingMode = fegetround();
   uint32_t r = ii_MaxN/Q - ii_MinN/Q + 1;
   
   SierpinskiRieselApp *srApp = (SierpinskiRieselApp *) ip_App;
   double   giantStepFactor = srApp->GetGiantStepFactor();
      
   fesetround(FE_TONEAREST);
   
   // r = range of n, divided by Q
   // s = number of subsequences.
   // In the worst case we will do do one table insertion and one mulmod
   // for m baby steps, then s table lookups and s mulmods for M giant
   // steps. The average case depends on how many solutions are found
   // and how early in the loop they are found, which I don't know how
   // to analyse. However for the worst case we just want to minimise
   // m + s*M subject to m*M >= r, which is when m = sqrt(s*r).

   M = MAX(1, rint(sqrt((double) (giantStepFactor * r)/s)));
   m = MIN(r, ceil((double) r/M));

   fesetround(roundingMode);
   
   // For either m or M to exceed this limit it would be really unusual.  It could only
   // happen if the range of N would have to be at least 2^30, s = 1 and Q = 1.
   if (m > SMALL_HASH_MAX_ELTS)
   {   
     // There are three ways to handle this undesirable case:
     //   1. Allow m > maxHashElements (undersize hash table).
     //   2. Restrict m <= maxHashElements (suboptimal baby/giant-step ratio).
     M = ceil((double)r/SMALL_HASH_MAX_ELTS);
     m = ceil((double)r/M);
   }
   
   while (m > SMALL_HASH_MAX_ELTS)
   {
      M++;
      m = ceil((double)r/M);
   }
   
   babySteps = m;
   giantSteps = M;
}