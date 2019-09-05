/* CIsOneSubsequenceHelper.cpp -- (C) Mark Rodenkirch, May 2019

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <assert.h>
#include "CIsOneSubsequenceHelper.h"
#include "NoLegendreWorker.h"
#include "LegendreWorker.h"
#include "HashTable.h"
#include "../core/inline.h"

#define NBIT(n)         ((n) - ii_MinN)
#define MBIT(m)         ((m) - ii_MinM)

CIsOneSubsequenceHelper::CIsOneSubsequenceHelper(App *theApp, seq_t *sequences, uint32_t sequenceCount, bool useLegendreTables, string legendreFileName) : AbstractSubsequenceHelper(theApp, sequences, sequenceCount)
{
   // The parent's constuctor does much of the work
   
   ib_UseLegendreTables = useLegendreTables;
   is_LegendreFileName = legendreFileName;

   theApp->WriteToConsole(COT_OTHER, "Sieving where abs(c) = 1 for all sequences");
}

Worker   *CIsOneSubsequenceHelper::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   NoLegendreWorker *theWorker;

   // Note that this inherits from Worker.  This will not
   // only create the worker, but also start it.
   if (ib_UseLegendreTables)
      theWorker = new LegendreWorker(id, ip_App, this);
   else
      theWorker = new NoLegendreWorker(id, ip_App, this);
   
   theWorker->SetSequences(largestPrimeTested, ii_BestQ, ip_Sequences, ii_SequenceCount, ip_Subsequences, ii_SubsequenceCount);
   
   return theWorker;
}

// Most of the logic below is copied from sr2sieve with changes to fit into 
// this application.  Most notably parameters and variable names have changed.

// Allow sieving in base b^Q for Q chosen from the divisors of LIMIT_BASE.
// LIMIT_BASE must be a multiple of BASE_MULTIPLE.

// These values are from sr2sieve.h
// For a prime p that satisfies p=1 (mod r), an "r-th power residue test"
// checks whether a subsequence of k*b^n+c can possibly contain any terms of
// the form x^r (mod p). If there are none then that subsequence can be
// omitted from the BSGS step.
//
// To conduct r-th power residue tests for each r in a set R of prime
// powers, set POWER_RESIDUE_LCM to lcm(R), and set POWER_RESIDUE_DIVISORS
// to the number of divisors of POWER_RESIDUE_LCM. POWER_RESIDUE_LCM must be
// a multiple of BASE_MULTIPLE, a divisor of LIMIT_BASE, and must be less
// than 2^15. E.g.
//
// R={2,3,4,5}: POWER_RESIDUE_LCM=60, POWER_RESIDUE_DIVISORS=12
// R={2,3,4,5,8}: POWER_RESIDUE_LCM=120, POWER_RESIDUE_DIVISORS=16
// R={2,3,4,5,8,9}: POWER_RESIDUE_LCM=360, POWER_RESIDUE_DIVISORS=24
// R={2,3,4,5,8,9,16}: POWER_RESIDUE_LCM=720, POWER_RESIDUE_DIVISORS=30
// R={2,3,4,5,7,8,9,16}: POWER_RESIDUE_LCM=5040, POWER_RESIDUE_DIVISORS=60
//
// Memory/sequence is proportional to POWER_RESIDUE_LCM*POWER_RESIDUE_DIVISORS

#define BASE_MULTIPLE      2
#define LIMIT_BASE         720
#define NDIVISORS          (LIMIT_BASE/BASE_MULTIPLE)
#define POWER_RESIDUE_LCM  720

void        CIsOneSubsequenceHelper::MakeSubsequencesForOldSieve(seq_t *sequences, uint64_t expectedTerms)
{
   uint32_t  j, r, n;
   uint32_t  bit, seqIdx, ssIdx;
   uint32_t  expectedSubsequences;
   uint32_t *rss;
   uint64_t  countedTerms = 0;
   bool     *needss;
   
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
   
   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      seq_t *seq = &ip_Sequences[seqIdx];

      for (j=0; j<ii_BestQ; j++)
         needss[j] = false;
      
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
   
   if (ii_SubsequenceCount > 1)
   {
      const char *sequenceText = ((ii_SequenceCount > 1) ? "sequences" : "sequence");
      
      ip_App->WriteToConsole(COT_OTHER, "Split %u base %u %s into %u base %u^%u sequences.", ii_SequenceCount,
            ii_Base, sequenceText, ii_SubsequenceCount, ii_Base, ii_BestQ);
   }
   
   for (ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      steps_t steps;
      
      ChooseSteps(ii_BestQ, ssIdx, steps);
      
      ip_Subsequences[ssIdx].babySteps = steps.babySteps;
      ip_Subsequences[ssIdx].giantSteps = steps.giantSteps;
   }
   
   MakeSubseqCongruenceTables();
      
   if (ib_UseLegendreTables)
   {
      if (is_LegendreFileName.length() > 0)
         LoadLegendreTables();
      else
         BuildLegendreTables();
   }
}

uint32_t    CIsOneSubsequenceHelper::FindBestQ(seq_t *sequences, uint32_t &expectedSubsequences)
{
   uint32_t i, j, n;
   uint32_t seqIdx, bit;
   uint32_t S[NDIVISORS];
   double   W[NDIVISORS];
   uint8_t  R[LIMIT_BASE];

   for (j=0; j<NDIVISORS; j++)
      S[j] = 0;

   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      seq_t *seq = &ip_Sequences[seqIdx];
      
      for (j=0; j<LIMIT_BASE; j++)
         R[j] = 0;
      
      bit = NBIT(ii_MinN);

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         if (seq->nTerms[bit])
            R[n%LIMIT_BASE] = 1;
         
         bit++;
      }
 
      for (j=0; j<NDIVISORS; j++)
      {
         if (NDIVISORS % (j+1) == 0)
            S[j] += CountResidueClasses((j+1)*BASE_MULTIPLE, LIMIT_BASE, R);
      }
   }

   for (i=0, j=0; j<NDIVISORS; j++)
   {
      if (NDIVISORS % (j+1) == 0)
      {
         W[j] = RateQ((j+1)*BASE_MULTIPLE, S[j]);
         if (W[j] < W[i])
            i = j;
      }
   }

   expectedSubsequences = S[i];
   
   return (i+1)*BASE_MULTIPLE;
}

uint32_t    CIsOneSubsequenceHelper::CountResidueClasses(uint32_t d, uint32_t Q, uint8_t *R)
{
   uint32_t i, count;
   uint8_t  R0[LIMIT_BASE];

   assert(Q % d == 0);

   for (i = 0; i < d; i++)
      R0[i] = 0;

   for (i = 0; i < Q; i++)
      if (R[i])
         R0[i%d] = 1;

   for (i = 0, count = 0; i < d; i++)
      if (R0[i])
         count++;

   return count;
}
                               
uint32_t    CIsOneSubsequenceHelper::RateQ(uint32_t Q, uint32_t s)
{
   double work, W[POWER_RESIDUE_LCM+1];
   uint32_t i;

   assert(Q % 2 == 0);
   assert(Q % BASE_MULTIPLE == 0);
   assert(LIMIT_BASE % Q == 0);

   for (i = 2, work = 0.0; i <= POWER_RESIDUE_LCM; i += 2)
   {
      if (POWER_RESIDUE_LCM % i == 0)
         W[i] = EstimateWork(Q, s, i);

      if (gcd32(i+1,POWER_RESIDUE_LCM) == 1)
         work += W[gcd32(i, POWER_RESIDUE_LCM)];
   }

   return work;
}

// These values are from choose.c
#define BABY_WORK    1.1    // 1 mulmod, 1 insert
#define GIANT_WORK   1.0    // 1 mulmod, 1 lookup
#define EXP_WORK     0.5    // 1 mulmod at most
#define SUBSEQ_WORK  1.0    // 1 mulmod, 1 lookup (giant step 0)
#define PRT_WORK     0.8    // List traversal, linear search, etc

// Return an estimate of the work needed to do one BSGS iteration on s
// subsequences in base b^Q with prime p%r=1.
double   CIsOneSubsequenceHelper::EstimateWork(uint32_t Q, uint32_t s, uint32_t r)
{
   steps_t steps;
   double work;
   uint32_t x;

   x = (s+r-1)/r; // subsequences expected to survive power residue tests
   
   ChooseSteps(Q, x, steps);
   
   work = steps.babySteps*BABY_WORK + x*(steps.giantSteps-1)*GIANT_WORK + Q*EXP_WORK + x*SUBSEQ_WORK;

   if (r > 2)
      work += x*PRT_WORK;

   return work;
}

void  CIsOneSubsequenceHelper::ChooseSteps(uint32_t Q, uint32_t s, steps_t &steps)
{
   uint32_t m, M;
   uint32_t r;

   r = ii_MaxN/Q - ii_MinN/Q+1;

   // r = range of n, s = number of subsequences.
   // In the worst case we will do do one table insertion and one mulmod
   // for m baby steps, then s table lookups and s mulmods for M giant
   // steps. The average case depends on how many solutions are found
   // and how early in the loop they are found, which I don't know how
   // to analyse. However for the worst case we just want to minimise
   // m + s*M subject to m*M >= r, which is when m = sqrt(s*r).

   M = MAX(1,rint(sqrt((double)r/s)));
   m = MIN(r,ceil((double)r/M));

   if (m > HASH_MAX_ELTS)
   {
      // There are three ways to handle this undesirable case:
      //   1. Allow m > HASH_MAX_ELTS (undersize hash table).
      //   2. Restrict m <= HASH_MAX_ELTS (suboptimal baby/giant-step ratio).
      //   3. Recompile with SHORT_HASHTABLE=0.
      M = ceil((double)r/HASH_MAX_ELTS);
      m = ceil((double)r/M);
   }

   assert(m <= HASH_MAX_ELTS);

   steps.babySteps = m;
   steps.giantSteps = M;
}

// Build tables sc_lists[i][j] of pointers to lists of subsequences
// whose terms k*b^m+c satisfy m = j (mod r) for i = divisor_index[r].
void  CIsOneSubsequenceHelper::MakeSubseqCongruenceTables(void)
{
   uint32_t  h, j, k, r, len, mem = 0;
   uint32_t  seqIdx, ssIdx;
   uint32_t  subseqList[POWER_RESIDUE_LCM];
   uint32_t *tmp;

   for (seqIdx=0; seqIdx<ii_SequenceCount; seqIdx++)
   {
      seq_t *seq = &ip_Sequences[seqIdx];
      
      for (j=1, r=0; j<=POWER_RESIDUE_LCM; j++)
      {
         if (POWER_RESIDUE_LCM % j == 0)
         {
            seq->sc_lists[r] = (uint32_t **) xmalloc(j*sizeof(uint32_t **));
            
            mem += j*sizeof(uint32_t *);
            
            for (k=0; k<j; k++)
            {
               for (len=0, ssIdx=seq->ssIdxFirst; ssIdx<=seq->ssIdxLast; ssIdx++)
                  if (CongruentTerms(ssIdx, j, k))
                  {
                     subseqList[len] = ssIdx;
                     len++;
                  }

               if (len == 0)
                  tmp = NULL;
               else
               {
                  tmp = (uint32_t *) xmalloc((len+1)*sizeof(uint32_t));
                  
                  mem += (len+1)*sizeof(uint32_t);
                  
                  for (h = 0; h < len; h++)
                     tmp[h] = subseqList[h];
                  
                  tmp[h] = UINT32_MAX;
               }
               
               seq->sc_lists[r][k] = tmp;
            }
            
            r++;
         }
      }
      
      assert(r == POWER_RESIDUE_DIVISORS);
   }
}

// Return true iff subsequence h of k*b^n+c has any terms with n%a==b.
bool  CIsOneSubsequenceHelper::CongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b)
{
   uint32_t g, m;
   subseq_t *subseq = &ip_Subsequences[ssIdx];

   g = gcd32(a, ii_BestQ);
   
   if (b % g == subseq->d % g)
   {
      for (m=ii_MinM; m<=ii_MaxM; m++)
      {        
         if (N_TERM(ssIdx, m) % a == b)
            return true;
      }
   }

   return false;
}