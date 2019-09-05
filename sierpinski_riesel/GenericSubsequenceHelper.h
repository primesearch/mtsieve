/* GenericSubsequenceHelper.h -- (C) Mark Rodenkirch, May 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) > 1.
*/

#ifndef _GenericSubsequenceHelper_H
#define _GenericSubsequenceHelper_H

#include "AbstractSubsequenceHelper.h"

class GenericSubsequenceHelper : public AbstractSubsequenceHelper
{
public:
   GenericSubsequenceHelper(App *theApp, seq_t *sequences, uint32_t sequenceCount);

   ~GenericSubsequenceHelper(void) {};

   void              MakeSubsequencesForOldSieve(seq_t *sequences, uint64_t expectedTermCount);

   uint64_t          SetSubsequenceTerms(seq_t *sequences, uint32_t sequenceCount);
   
   void              ComputeSmallSieveLimit(void);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
private:
   uint32_t          FindBestQ(seq_t *sequences, uint32_t &expectedSubsequences);
   
   uint32_t          ForEachDivisor(uint32_t Q, bool *R, choice_bc_t *S);
   uint32_t          FindMultiplicities(uint32_t n, uint32_t *P, uint32_t *M);
   uint32_t          CountResidueClasses(uint32_t d, uint32_t Q, bool *R);
   uint32_t          RateQ(uint32_t Q, uint32_t s);
   void              ChooseSteps(uint32_t *baby, uint32_t *giant, uint32_t Q, uint32_t s);
};

#endif

