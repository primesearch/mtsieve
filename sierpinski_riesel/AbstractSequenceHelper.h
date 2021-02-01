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

// Set BASE_MULTIPLE to the smallest exponent Q for which sieving in
// subsequence base b^Q will be considered. Must be a multiple of 2.
#define BASE_MULTIPLE   30

// Allow sieving in base b^Q for Q chosen from the divisors of LIMIT_BASE.
// Must be a multiple of BASE_MULTIPLE.
#define LIMIT_BASE         720

#define NDIVISORS          (LIMIT_BASE/BASE_MULTIPLE)

// For a prime p that satisfies p=1 (mod r), an "r-th power residue test"
// checks whether a subsequence of k*b^n+c can possibly contain any terms of
// the form x^r (mod p). If there are none then that subsequence can be
// omitted from the BSGS step.
//
// To conduct r-th power residue tests for each r in a set R of prime powers,
// set POWER_RESIDUE_LCM to lcm(R). POWER_RESIDUE_LCM must be a multiple of
// BASE_MULTIPLE, a divisor of LIMIT_BASE, and must be less than 2^15.
#define POWER_RESIDUE_LCM      720

#define SP_COUNT   3
typedef enum { SP_NO_PARITY = 999, SP_MIXED = 0, SP_EVEN = 1, SP_ODD = 2} sp_t;

// This allows us to create a one dimensional array to access the qList and ladders.
// Note that this goes from the largest dimension to the smallest so that each x/y/z
// combination yeiel
#define CSS_INDEX(x, y, z) (((((x) * ii_PrlCount) + (y)) * POWER_RESIDUE_LCM) + (z))

// The maps are not vector<bool> because if I ever write an OpenCL
// kernel for this, it has to be a simple datatype.
typedef struct {
   uint8_t          *oneParityMap;
   uint8_t          *dualParityMapM1;
   uint8_t          *dualParityMapP1;
   uint32_t          mod;
} legendre_t;

// All of these fields are set before sieving is started, but only nTerms can be
// modified after sieving has started.
typedef struct
{
   uint64_t     k;            // k in (k*b^n+c)/d
   int64_t      c;            // c in (k*b^n+c)/d
   uint32_t     d;            // d in (k*b^n+c)/d
   
   uint64_t     squareFreeK;  // used by the CIsOne classes, product of square free factors of k
   int64_t      kcCore;       // used by the CIsOne classes, squareFreeK * -c, only used if unable to use Legendre tables
   
   sp_t         parity;       // 0 if all n terms are even, 1 if all n terms are odd, 2 if mixed
   uint32_t     ssCount;
   uint32_t     ssIdxFirst;   // index of first subsequence for the sequence
   uint32_t     ssIdxLast;    // index of first subsequence for the sequence
   
   vector<bool> nTerms;       // remaining n for this sequences
   
   MpResVec     resCK;        // scratch space used by the workers
   void        *next;         // points to the next sequence

   // The fields below are only used by the CisOne classes
   legendre_t  *legendrePtr; 
   
   // Congruent subsequence details
   // In sr1sieve, these handled via four dimensional arrays.
   // For srsieve2, we will use two one dimensional arrays, which will be easier to pass to the GPU.
   // congruentQIndices points to the first entry in congruentQs for a specific parity, r, and h.
   // congruentQs is a list of qs for that parity, r, and h with the first entry the length of the list
   // for that pariry, r, and h.  The relationship for the ladder is the same.
   // This will also reduce the memory needed by the GPU to hold these structures.
   
   uint32_t    *congruentQIndices;
   uint32_t    *congruentLadderIndices;
   
   uint32_t     congruentQSize;
   uint32_t     congruentLadderSize;
   
   uint16_t    *congruentQs;
   uint16_t    *congruentLadders;
} seq_t;

// All of these fields are set before sieving is started, but none of them are
// modified after sieving has started.
typedef struct
{
   seq_t       *seqPtr;       // points to the seq_t that this subsequence is for
   uint64_t     k;            // k in k*b^n+c
   int64_t      c;            // c in k*b^n+c
   uint32_t     q;
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

   // These are only used by the CisOne logic
   uint32_t          GetPrlCount(void) { return ii_PrlCount; };
   uint16_t         *GetPrlIndices(void) { return ip_PrlIndices; };
   
protected:
   void              CreateEmptySubsequences(uint32_t subsequenceCount);
   uint32_t          AddSubsequence(seq_t *seq, uint32_t q, uint32_t mTermCount);
   
   virtual uint32_t  FindBestQ(uint32_t &expectedSubsequences) = 0;
   
   uint32_t          CountResidueClasses(uint32_t d, uint32_t Q, bool *R);
   void              ChooseSteps(uint32_t Q, uint32_t s, uint32_t &babySteps, uint32_t &giantSteps);
   
   App              *ip_App;
   
   seq_t            *ip_FirstSequence;
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

   // These are only used by the CisOne logic.
   // Only a small subset of values < POWER_RESIDUE_LCM can generate values to populate
   // the qList and ladders.  This is used to reduce the amount of memory needed for
   // seq->congruentQIndices and seq->congruentLadderIndices, which will benefit the GPU.
   uint32_t          ii_PrlCount;
   uint16_t         *ip_PrlIndices;

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

