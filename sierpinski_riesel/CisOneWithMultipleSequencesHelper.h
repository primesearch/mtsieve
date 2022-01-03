/* CisOneWithMultipleSequencesHelper.h -- (C) Mark Rodenkirch, January 2021
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This class is used if any sequence has abs(c) = 1.
*/

#ifndef _CisOneWithMultipleSequencesHelper_H
#define _CisOneWithMultipleSequencesHelper_H

#include "CisOneSequenceHelper.h"

// This allows us to create a one dimensional array to access the qList and ladders
//    a = seqIdx
//    b = power residue index
//    c = power residue for congruent q
#define CSS_INDEX(a, b, c) ((a) * ii_Dim2 + (b) * ii_Dim3 + (c))

#define DEFAULT_BM_MULTIPLIER_MULTI    1
#define DEFAULT_PRL_MULTIPLIER_MULTI   360
#define DEFAULT_LB_MULTIPLIER_MULTI    1

class CisOneWithMultipleSequencesHelper : public CisOneSequenceHelper
{
public:
   CisOneWithMultipleSequencesHelper(App *theApp, uint64_t largestPrimeTested);

   ~CisOneWithMultipleSequencesHelper(void) {};
   
   Worker        *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   uint32_t      *GetCongruentSubseqIndices(void) { return ip_CongruentSubseqIndices; };
   uint32_t      *GetAllSubseqs(void) { return ip_AllSubseqs; };
   uint16_t      *GetAllLadders(void) { return ip_AllLadders; };

   uint32_t       GetDim1(void) { return ii_Dim1; };
   uint32_t       GetDim2(void) { return ii_Dim2; };
   uint32_t       GetDim3(void) { return ii_Dim3; };
   
protected:
   double         RateQ(uint32_t Q, uint32_t s);
   double         EstimateWork(uint32_t Q, uint32_t s, uint32_t r);

   void           BuildCongruenceTables(void);

   void           CopySubseqs(seq_t *seqPtr, uint32_t r, uint32_t h, uint32_t *ssList, uint32_t ssListLen);
   void           MakeLadder(void);

   void           ComputeLegendreMemoryToAllocate(legendre_t *legendrePtr, uint64_t ssqfb);
   void           AssignMemoryToLegendreTable(legendre_t *legendrePtr, uint64_t bytesUsed);
   void           BuildLegendreTableForSequence(legendre_t *legendrePtr, uint64_t ssqfb, uint64_t stepsToDo, uint64_t stepsDone, time_t startTime);

   uint32_t       ii_Dim1;
   uint32_t       ii_Dim2;
   uint32_t       ii_Dim3;
   
   uint32_t       ii_UsedSubseqEntries;
   uint32_t       ii_LadderEntries;
   uint32_t       ii_MaxSubseqEntries;
   
   uint32_t      *ip_CongruentSubseqIndices;
   
   uint32_t      *ip_AllSubseqs;
   uint16_t      *ip_AllLadders;
};

#endif

