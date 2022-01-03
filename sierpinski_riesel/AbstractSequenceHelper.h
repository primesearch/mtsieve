/* AbstractSequenceHelper.h -- (C) Mark Rodenkirch, May 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AbstractSubsequenceHelper_H
#define _AbstractSubsequenceHelper_H

#include "../core/Worker.h"
#include "../core/MpArithVector.h"

#define SP_COUNT   3
typedef enum { SP_NO_PARITY = 999, SP_MIXED = 0, SP_EVEN = 1, SP_ODD = 2} sp_t;

// All of these fields are set before sieving is started, but only nTerms can be
// modified after sieving has started.
typedef struct
{
   uint32_t     seqIdx;       // sequence index, not used for the global sequences

   uint64_t     k;            // k in (k*b^n+c)/d
   int64_t      c;            // c in (k*b^n+c)/d
   uint32_t     d;            // d in (k*b^n+c)/d
   
   uint64_t     squareFreeK;  // used by the CIsOne classes, product of square free factors of k
   int64_t      kcCore;       // used by the CIsOne classes, squareFreeK * -c, only used if unable to use Legendre tables
   
   sp_t         nParity;      // 0 if all n terms are even, 1 if all n terms are odd, 2 if mixed
   uint32_t     ssCount;
   uint32_t     ssIdxFirst;   // index of first subsequence for the sequence
   uint32_t     ssIdxLast;    // index of first subsequence for the sequence
   
   vector<bool> nTerms;       // remaining n for this sequences
   
   void        *next;         // points to the next sequence
} seq_t;

// All of these fields are set before sieving is started, but none of them are
// modified after sieving has started.
typedef struct
{
   seq_t       *seqPtr;       // points to the seq_t that this subsequence is for
   uint64_t     k;            // k in k*b^n+c
   int64_t      c;            // c in k*b^n+c
   uint16_t     q;
   vector<bool> mTerms;       // remaining m for this sub-sequence
   
   uint32_t     babySteps;    // baby steps for CIsOne logic
   uint32_t     giantSteps;   // giant steps for CIsOne logic
} subseq_t;

typedef struct
{
  uint32_t div;
  uint32_t subseqs;
  double   work;
} choice_bc_t;

class AbstractSequenceHelper
{
public:
   AbstractSequenceHelper(App *theApp, uint64_t largestPrimeTested);
   
   virtual ~AbstractSequenceHelper(void) {};

   void              CleanUp(void);
      
   uint64_t          MakeSubsequencesForNewSieve(void);
   
   void              MakeSubsequencesForOldSieve(uint64_t expectedTermCount);

   virtual void      LastChanceLogicBeforeSieving(void) = 0;
   
   virtual Worker   *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested) = 0;
   
   void              NotifyWorkerOfUpdatedSequences(Worker *worker);

   seq_t            *GetFirstSequenceAndSequenceCount(uint32_t &count) { count = ii_SequenceCount; return ip_FirstSequence; };
   subseq_t         *GetSubsequences(uint32_t &count) { count = ii_SubsequenceCount; return ip_Subsequences; };
   
protected:
   void              CreateEmptySubsequences(uint32_t subsequenceCount);
   uint32_t          AddSubsequence(seq_t *seqPtr, uint32_t q, uint32_t mTermCount);
   
   virtual uint32_t  FindBestQ(uint32_t &expectedSubsequences) = 0;
   
   uint32_t          CountResidueClasses(uint32_t d, uint32_t Q, vector<bool> R);
   void              ChooseSteps(uint32_t Q, uint32_t s, uint32_t &babySteps, uint32_t &giantSteps);
   
   App              *ip_App;
   
   seq_t            *ip_FirstSequence;
   subseq_t         *ip_Subsequences;
   
   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint32_t          ii_SequenceCount;
   uint32_t          ii_MaxSubsequenceCount;
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

