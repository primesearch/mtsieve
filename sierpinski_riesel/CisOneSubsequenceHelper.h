/* CIsOneSubsequenceHelper.h -- (C) Mark Rodenkirch, May 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) = 1.
*/

#ifndef _CIsOneSubsequenceHelper_H
#define _CIsOneSubsequenceHelper_H

#include "AbstractSubsequenceHelper.h"

typedef struct {
   uint32_t babySteps; 
   uint32_t giantSteps;
} steps_t;

class CIsOneSubsequenceHelper : public AbstractSubsequenceHelper
{
public:
   CIsOneSubsequenceHelper(App *theApp, seq_t *sequences, uint32_t sequenceCount, string legendreFileName);

   ~CIsOneSubsequenceHelper(void) {};

   void              MakeSubsequencesForOldSieve(seq_t *sequences, uint64_t expectedTermCount);

   uint64_t          SetSubsequenceTerms(seq_t *sequences, uint32_t sequenceCount);
   
   void              ComputeSmallSieveLimit(void);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
private:
   bool              ib_UseLegendreTables;
   string            is_LegendreFileName;
   
   uint32_t          FindBestQ(seq_t *sequences, uint32_t &expectedSubsequences);
   
   uint32_t          ForEachDivisor(uint32_t Q, bool *R, choice_bc_t *S);
   uint32_t          FindMultiplicities(uint32_t n, uint32_t *P, uint32_t *M);
   uint32_t          CountResidueClasses(uint32_t d, uint32_t Q, uint8_t *R);
   uint32_t          RateQ(uint32_t Q, uint32_t s);
   void              ChooseSteps(uint32_t Q, uint32_t s, steps_t &steps);
   void              MakeSubseqCongruenceTables(void);
   double            EstimateWork(uint32_t Q, uint32_t s, uint32_t r);
   bool              CongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b);
};

#endif

