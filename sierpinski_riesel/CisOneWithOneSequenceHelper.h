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

// This allows us to create a one dimensional array to access the qList and ladders
//    a = parity
//    b = power residue index
//    c = power residue for congruent q
#define CQ_INDEX(a, b, c) ((a) * ii_Dim2 + (b) * ii_Dim3 + (c))

#define DEFAULT_BM_MULTIPLIER_SINGLE   15
#define DEFAULT_PRL_MULTIPLIER_SINGLE  24
#define DEFAULT_LB_MULTIPLIER_SINGLE   1

class CisOneWithOneSequenceHelper : public CisOneSequenceHelper
{
public:
   CisOneWithOneSequenceHelper(App *theApp, uint64_t largestPrimeTested);

   ~CisOneWithOneSequenceHelper(void) {};
   
   Worker        *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

   uint32_t      *GetCongruentQIndices(void) { return ip_CongruentQIndices; };
   uint32_t      *GetLadderIndices(void) { return ip_LadderIndices; };
   uint16_t      *GetAllQs(void) { return ip_AllQs; };
   uint16_t      *GetAllLadders(void) { return ip_AllLadders; };

   uint32_t       GetUsedQEntries(void) { return ii_UsedQEntries; };
   uint32_t       GetUsedLadderEntries(void) { return ii_UsedLadderEntries; };

   uint32_t       GetDim1(void) { return ii_Dim1; };
   uint32_t       GetDim2(void) { return ii_Dim2; };
   uint32_t       GetDim3(void) { return ii_Dim3; };
   
protected:
   double         RateQ(uint32_t Q, uint32_t s);
   double         EstimateWork(uint32_t Q, uint32_t s);
   
   void           BuildCongruenceTables(void);
   void           BuildCongruenceTablesForSequence(seq_t *seqPtr);

   void           CopyQsAndMakeLadder(seq_t *seqPtr, sp_t parity, uint32_t r, uint32_t h, uint16_t *qList, uint32_t qListLen);
   void           MakeLadder(uint16_t *qList, uint32_t qListLen);
   
   void           ComputeLegendreMemoryToAllocate(legendre_t *legendrePtr, uint64_t ssqfb);
   void           AssignMemoryToLegendreTable(legendre_t *legendrePtr, uint64_t bytesUsed);
   void           BuildLegendreTableForSequence(legendre_t *legendrePtr, uint64_t ssqfb, uint64_t stepsToDo, uint64_t stepsDone, time_t startTime);

   uint32_t       ii_Dim1;
   uint32_t       ii_Dim2;
   uint32_t       ii_Dim3;
   
   uint32_t       ii_UsedQEntries;
   uint32_t       ii_UsedLadderEntries;
   uint32_t       ii_MaxQEntries;
   uint32_t       ii_MaxLadderEntries;
   
   uint32_t      *ip_CongruentQIndices;
   uint32_t      *ip_LadderIndices;
   
   uint16_t      *ip_AllQs;
   uint16_t      *ip_AllLadders;
};

#endif

