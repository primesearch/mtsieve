/* GenericSequenceHelper.h -- (C) Mark Rodenkirch, May 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) > 1.
*/

#ifndef _GenericSequenceHelper_H
#define _GenericSequenceHelper_H

#include "AbstractSequenceHelper.h"

class GenericSequenceHelper : public AbstractSequenceHelper
{
public:
   GenericSequenceHelper(App *theApp, uint64_t largestPrimeTested);

   ~GenericSequenceHelper(void) {};

   void              LastChanceLogicBeforeSieving(void) {};

   void              ComputeSmallSieveLimit(void);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
private:
   uint32_t          FindBestQ(uint32_t &expectedSubsequences);
   
   uint32_t          ForEachDivisor(uint32_t Q, std::vector<bool> R, choice_bc_t *S, bool firstTime);
   uint32_t          FindMultiplicities(uint32_t n, uint32_t *P, uint32_t *M);
   double            EstimateWork(uint32_t Q, uint32_t s);
};

#endif

