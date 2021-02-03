/* CisOneSequenceHelper.h -- (C) Mark Rodenkirch, May 2019
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) = 1.
*/

#ifndef _CisOneSequenceHelper_H
#define _CisOneSequenceHelper_H

#include "AbstractSequenceHelper.h"

class CisOneSequenceHelper : public AbstractSequenceHelper
{
public:
   CisOneSequenceHelper(App *theApp, uint64_t largestPrimeTested);

   ~CisOneSequenceHelper(void) {};

   void              CleanUp(void);
   
   void              LastChanceLogicBeforeSieving(void);
   
   void              ComputeSmallSieveLimit(void);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested) = 0;

   uint32_t          GetMaxBabySteps(void) { return ii_MaxBabySteps; };
   
   int16_t          *GetDivisorShifts(void) { return ip_DivisorShifts; };
   
   bool              BuildLegendreTables(string legendreFileName);
   bool              UseLegendreTables(void) { return ib_UseLegendreTables; };

protected:
   virtual bool      BuildLegendreTableForSequence(seq_t *seq) = 0;
   
   bool              ib_UseLegendreTables;
   string            is_LegendreFileName;

   uint32_t          FindBestQ(uint32_t &expectedSubsequences);
   virtual uint32_t  RateQ(uint32_t Q, uint32_t s) = 0;
   double            EstimateWork(uint32_t Q, uint32_t s, uint32_t r);
   
   virtual void      MakeSubseqCongruenceTables(seq_t *seq) = 0;

   uint32_t          ii_MaxBabySteps;
   uint32_t          ii_LadderCount;
   uint32_t          ii_MaxQs;
   uint32_t          ii_MaxLadderRungs;
      
   uint16_t         *ip_TempQs;
   uint16_t         *ip_TempLadders;
   
   uint32_t          ii_TempQIndex;
   uint32_t          ii_TempLadderIndex;
   
   int16_t          *ip_DivisorShifts;
};

#endif

