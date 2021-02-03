/* CisOneWithOneSequenceHelper.h -- (C) Mark Rodenkirch, January 2021
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) = 1.
*/

#ifndef _CisOneWithOneSequenceHelper_H
#define _CisOneWithOneSequenceHelper_H

#include "CisOneSequenceHelper.h"

class CisOneWithOneSequenceHelper : public CisOneSequenceHelper
{
public:
   CisOneWithOneSequenceHelper(App *theApp, uint64_t largestPrimeTested);

   ~CisOneWithOneSequenceHelper(void) {};
   
   Worker     *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
protected:
   uint32_t    RateQ(uint32_t Q, uint32_t s);
   double      EstimateWork(uint32_t Q, uint32_t s);

   bool        BuildLegendreTableForSequence(seq_t *seq);
   void        BuildLegendreMap(uint32_t size, int64_t r, const char *which, uint8_t *legendreMap);
   
   void        MakeSubseqCongruenceTables(seq_t *seq);
   bool        CongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b);

   bool        CopyQsAndMakeLadder(seq_t *seq, sp_t parity, uint32_t r, uint32_t h, uint16_t *qList, uint32_t qListLen);

   void        MakeLadder(seq_t *seq, uint16_t *qList, uint32_t qListLen);
};

#endif
