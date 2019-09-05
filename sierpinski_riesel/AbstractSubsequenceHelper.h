/* AbstractSubsequenceHelper.h -- (C) Mark Rodenkirch, May 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AbstractSubsequenceHelper_H
#define _AbstractSubsequenceHelper_H

#include "../core/Worker.h"

// from sr2sieve.h
#define POWER_RESIDUE_DIVISORS 30

#define SEQ_IDX(ssIdx)     (ip_Subsequences[(ssIdx)].seqIdx)
#define N_TERM(ssIdx, m)   ((m)*ii_BestQ + ip_Subsequences[(ssIdx)].d)

// All of these fields are set before sieving is started, but only nTerms can be
// modified after sieving has started.
typedef struct
{
   uint64_t     k;            // k in k*b^n+c
   int64_t      c;            // c in k*b^n+c
   
   uint32_t     parity;       // 0 if all n terms are even, 1 if all n terms are odd, 2 if mixed
   uint32_t     ssCount;
   uint32_t     ssIdxFirst;   // index of first subsequence for the sequence
   uint32_t     ssIdxLast;    // index of first subsequence for the sequence
   
   uint32_t   **sc_lists[POWER_RESIDUE_DIVISORS];
   
   vector<bool> nTerms;       // remaining n for this sequence
   
   uint64_t     legendreMod;  // used only if ib_UseLengendreTables is true
   vector<bool> legendreMap;  // used only if ib_UseLengendreTables is true
} seq_t;

// All of these fields are set before sieving is started, but none of them are
// modified after sieving has started.
typedef struct
{
   uint32_t     seqIdx;       // index to the sequence that this subsequence is for
   uint64_t     k;            // k in k*b^n+c
   int64_t      c;            // c in k*b^n+c
   uint32_t     d;
   uint32_t     babySteps;
   uint32_t     giantSteps;
   vector<bool> mTerms;       // remaining m for this sub-sequence
} subseq_t;

typedef struct
{
  uint32_t div;
  uint32_t subseqs;
  uint32_t work;
} choice_bc_t;

class AbstractSubsequenceHelper
{
public:
   AbstractSubsequenceHelper(App *theApp, seq_t *sequences, uint32_t sequenceCount);
   
   virtual ~AbstractSubsequenceHelper(void) {};

   void              CleanUp(void);
   
   uint32_t          GetSubsequenceCount(void) { return ii_SubsequenceCount; };
   
   uint64_t          MakeSubsequencesForNewSieve(seq_t *sequences);
   
   virtual void      MakeSubsequencesForOldSieve(seq_t *sequences, uint64_t expectedTermCount) = 0;

   virtual Worker   *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested) = 0;
   
   void              NotifyWorkerOfUpdatedSequences(Worker *worker);
      
protected:
   uint32_t          AddSubsequence(uint32_t seqIdx, uint32_t d, uint32_t mTermCount);
   
   App              *ip_App;
   
   seq_t            *ip_Sequences;
   subseq_t         *ip_Subsequences;
   
   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_SequenceCount;
   uint32_t          ii_SubsequenceCount;
   uint32_t          ii_SubsequenceCapacity;

   uint32_t          ii_BestQ;
   uint32_t          ii_MinM;
   uint32_t          ii_MaxM;
   
   inline uint32_t   pow32(uint32_t b, uint32_t n)
   {
      uint32_t a = 1;

      while (n > 0)
      {
         if (n % 2)
            a *= b;
         b *= b;
         n /= 2;
      }

      return a;
   }
};

#endif

