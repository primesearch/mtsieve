/* GenericSubsequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "GenericSubsequenceHelper.h"
#include "GenericWorker.h"

#ifdef HAVE_GPU_WORKERS
#include "GenericGpuWorker.h"
#endif

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

GenericSubsequenceHelper::GenericSubsequenceHelper(App *theApp, seq_t *sequences, uint32_t sequenceCount) : AbstractSubsequenceHelper(theApp, sequences, sequenceCount)
{
   // The parent's constuctor does the work
   theApp->WriteToConsole(COT_OTHER, "Sieving with generic logic");
}

Worker   *GenericSubsequenceHelper::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   AbstractWorker *theWorker;

   // Note that GenericWorker inherits from Worker.  This will not
   // only create the worker, but also start it.
   
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      theWorker = new GenericGpuWorker(id, ip_App, this);
   else
#endif
      theWorker = new GenericWorker(id, ip_App, this);
      
   theWorker->SetSequences(largestPrimeTested, ii_BestQ, ip_Sequences, ii_SequenceCount, ip_Subsequences, ii_SubsequenceCount);

   return theWorker;
}

void        GenericSubsequenceHelper::MakeSubsequencesForOldSieve(seq_t *sequences, uint64_t expectedTerms)
{
   bool     *needss;
   uint32_t  bit, r, n;
   uint32_t *rss;
   uint32_t  expectedSubsequences;
   uint64_t  countedTerms = 0;
   uint32_t  seqIdx, ssIdx;
   
   if (ip_Subsequences)
      xfree(ip_Subsequences);
   
   ip_Subsequences = 0;
   ii_SubsequenceCount = 0;
   ii_SubsequenceCapacity = 0;
   
   ii_BestQ = FindBestQ(sequences, expectedSubsequences);
   
   ii_MinM = ii_MinN/ii_BestQ;
   ii_MaxM = ii_MaxN/ii_BestQ;

   needss = (bool *) xmalloc(ii_BestQ*sizeof(bool));
   rss = (uint32_t *) xmalloc(ii_BestQ*sizeof(uint32_t));

   // This is the maximum number of subsequences
   CreateEmptySubsequences(ii_SequenceCount * ii_BestQ);
   
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      seq_t *seq = &ip_Sequences[seqIdx];
      
      for (r=0; r<ii_BestQ; r++)
         needss[r] = false;

      bit = NBIT(ii_MinN);
      
      // Determine which subsequences we need to build
      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[bit])
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
            AddSubsequence(seqIdx, r, ii_MaxM - ii_MinM + 1);
         }
      }
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[bit])
         {
            r = n % ii_BestQ;
            ssIdx = rss[r];
            ip_Subsequences[ssIdx].mTerms[n/ii_BestQ - ii_MinM] = true;
            countedTerms++;
         }
         
         bit++;
      }
   }

   xfree(rss);
   xfree(needss);

   if (ii_SubsequenceCount != expectedSubsequences)
      FatalError("Expected %u subsequences but %u were created", expectedSubsequences, ii_SubsequenceCount);

   if (expectedTerms != countedTerms)
      FatalError("Expected %" PRIu64" terms when building sequences, but counted only %" PRIu64"", expectedTerms, countedTerms);
   
  const char *sequenceText = ((ii_SequenceCount > 1) ? "sequences" : "sequence");
  
  ip_App->WriteToConsole(COT_OTHER, "Split %u base %u %s into %u base %u^%u sequences.", ii_SequenceCount,
		ii_Base, sequenceText, ii_SubsequenceCount, ii_Base, ii_BestQ);
}

uint32_t    GenericSubsequenceHelper::FindBestQ(seq_t *sequences, uint32_t &expectedSubsequences)
{
   uint32_t bit, i, j, k, n;
   uint32_t seqIdx;
   uint32_t bestQ, Q = 2880;  // DEFAULT_LIMIT_BASE from the old code
   bool    *R = 0;
   choice_bc_t *S = 0;

   k = ForEachDivisor(Q, R, S);
   
   R = (bool *) xmalloc(Q*sizeof(bool));
   S = (choice_bc_t *) xmalloc(k*sizeof(choice_bc_t));

   for (i = 0; i < k; i++)
      S[i].subseqs = 0;

   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      seq_t *seq = &ip_Sequences[seqIdx];
      
      for (j=0; j<Q; j++)
         R[j] = 0;
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[bit])
            R[n%Q] = 1;
         
         bit++;
      }
      
      ForEachDivisor(Q, R, S);
   }
   
   for (i = 0; i < k; i++)
      S[i].work = RateQ(S[i].div, S[i].subseqs);

   j = 0;
   for (i = 0; i < k; i++)
      if (S[i].work < S[j].work)
         j = i;
   
   bestQ = S[j].div;
   expectedSubsequences = S[j].subseqs;
   
   xfree(R);
   xfree(S);
   
   return bestQ;
}

uint32_t GenericSubsequenceHelper::ForEachDivisor(uint32_t Q, bool *R, choice_bc_t *S)
{
   uint32_t A[9], P[9], M[9], d, i, j, k, t;

   k = FindMultiplicities(Q, P, M);

   for (i = 0, t = 1; i < k; i++)
   {
      A[i] = 0;
      t *= (M[i]+1);
   }

   if (!R)
      return t;
   
   for (i = 0; i < t; i++)
   {
      for (j = 0; j < k; A[j++] = 0)
         if (++A[j] <= M[j])
            break;
         
      for (j = 0, d = 1; j < k; j++)
         d *= pow32(P[j], A[j]);

     S[i].div = d;
     S[i].subseqs += CountResidueClasses(d, Q, R);
   }

   return t;
}

// Express n as n = p_0^m_0 * p_1^m_1 * ... * p_(k-1)^m_(k-1) for prime p_i
// and positive m_j. Store p_i in P[i], m_j in M[j], and return k. Since n
// is a uint32_t, k cannot be greater than 9.
uint32_t    GenericSubsequenceHelper::FindMultiplicities(uint32_t n, uint32_t *P, uint32_t *M)
{
   uint32_t i, m, q;

   for (i = 0, q = 2; n > 1; q = (q + 1) | 1)
   {
      for (m = 0; n % q == 0; m++)
         n /= q;
      
      if (m > 0)
         P[i] = q, M[i] = m, i++;
   }

   return i;
}

uint32_t    GenericSubsequenceHelper::CountResidueClasses(uint32_t d, uint32_t Q, bool *R)
{
   uint32_t i, count = 0;
   bool       *R0;

   R0 = (bool *) xmalloc(d*sizeof(bool));

   for (i = 0; i < d; i++)
      R0[i] = 0;

   for (i = 0; i < Q; i++)
      if (R[i])
         R0[i%d] = 1;

   for (i = 0; i < d; i++)
      if (R0[i])
         count++;

   xfree(R0);

   return count;
}

#define BABY_WORK    1.0    // 1 mulmod, 1 insert
#define GIANT_WORK   1.0    // 1 mulmod, 1 lookup
#define EXP_WORK     0.7    // 1 mulmod
#define SUBSEQ_WORK  1.4    // 1 mulmod, 1 lookup (giant step 0)
                               
uint32_t    GenericSubsequenceHelper::RateQ(uint32_t Q, uint32_t s)
{
   uint32_t baby, giant, work;

   ChooseSteps(&baby, &giant, Q, s);

   work = baby*BABY_WORK + s*(giant-1)*GIANT_WORK + Q*EXP_WORK + s*SUBSEQ_WORK;

   return work;
}

void     GenericSubsequenceHelper::ChooseSteps(uint32_t *baby, uint32_t *giant, uint32_t Q, uint32_t s)
{
   uint32_t r = ii_MaxN/Q - ii_MinN/Q + 1;

   *giant = MAX(1,sqrt((double)r/s));
   *baby = MIN(r,ceil((double)r/(*giant)));
}
