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

// The maps are not vector<bool> because if I ever write an OpenCL
// kernel for this, it has to be a simple datatype.
typedef struct {
   uint8_t          *oneParityMap;
   uint8_t          *dualParityMapM1;
   uint8_t          *dualParityMapP1;
   uint32_t          mod;
} legendre_t;

class CisOneSequenceHelper : public AbstractSequenceHelper
{
public:
   CisOneSequenceHelper(App *theApp, uint64_t largestPrimeTested);

   ~CisOneSequenceHelper(void) {};

   void              CleanUp(void);
   
   void              LastChanceLogicBeforeSieving(void);
   
   void              ComputeSmallSieveLimit(void);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested) = 0;

   uint32_t          GetSieveLow(void) { return ii_SieveLow; };
   uint32_t          GetSieveHigh(void) { return ii_SieveHigh; };
   uint32_t          GetMaxBabySteps(void) { return ii_MaxBabySteps; };
   
   int16_t          *GetDivisorShifts(void) { return ii_DivisorShifts; };
   uint16_t         *GetCongruenceList(uint32_t parity, uint32_t r, uint32_t h) { return ii_CssQs[parity][r][h]; };
   uint16_t         *GetCongruenceLadder(uint32_t parity, uint32_t r, uint32_t h) { return ii_CSSLadders[parity][r][h]; };

   bool              BuildLegendreTables(string legendreFileName);
   bool              UseLegendreTables(void) { return ib_UseLegendreTables; };

protected:
   virtual bool      BuildLegendreTableForSequence(seq_t *seq, legendre_t *legendre) = 0;
   
   bool              ib_UseLegendreTables;
   string            is_LegendreFileName;

   uint32_t          FindBestQ(uint32_t &expectedSubsequences);
   virtual uint32_t  RateQ(uint32_t Q, uint32_t s) = 0;
   double            EstimateWork(uint32_t Q, uint32_t s, uint32_t r);
   
   void              MakeSubseqCongruenceTables(void);
   bool              CongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b);
   uint16_t         *MakeLadder(uint16_t *ssList, uint32_t len);
   
   uint32_t          ii_SieveLow;
   uint32_t          ii_SieveHigh;
   uint32_t          ii_MaxBabySteps;
   
   uint64_t          ii_MaxQs;
   uint32_t          ii_MaxLadderRungs;
   int16_t          *ii_DivisorShifts;
   
   // Congruent subsequence details
   uint16_t       ***ii_CssQs[3];
   uint16_t       ***ii_CSSLadders[3];
};

#endif

