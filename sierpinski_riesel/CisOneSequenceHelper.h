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

#define L_BYTES(x) (((1+x)>>3)+1)
#define L_BYTE(x)  ((x)>>3)
#define L_BIT(x)   (1<<((x)&7))

// The maps are not vector<bool> because it has to be a simple datatype for OpenCL
typedef struct {
   uint32_t    seqIdx;
   uint64_t    k;
   int64_t     c;
   int64_t     kcCore;
   sp_t        nParity;
   uint64_t    squareFreeK;
   uint32_t    mod;
   int64_t     r;
   
   uint32_t    bytesNeeded;            // total space needed for all maps
   bool        canCreateMap;           // indicates if we can build Ledendre tables for this sequence
   bool        haveMap;                // indicates if we have built the Legendre tables for this sequence
   bool        loadedMapFromCache;
   
   uint32_t    mapSize;                // size of each map in bytes
   
   uint8_t    *oneParityMap;
   uint8_t    *dualParityMapM1;        // not used for CisOneWithMultiSequenceHelper
   uint8_t    *dualParityMapP1;        // not used for CisOneWithMultiSequenceHelper
  
   uint64_t    oneParityMapIndex;      // used by the GPU
   uint64_t    dualParityMapM1Index;   // used by the GPU, not used for CisOneWithMultiSequenceHelper
   uint64_t    dualParityMapP1Index;   // used by the GPU, not used for CisOneWithMultiSequenceHelper
} legendre_t;

typedef struct {
   uint32_t    fileVersion;
   uint32_t    structSize;

   uint32_t    base;
   uint64_t    k;
   int64_t     c;
   
   uint32_t    mapSize;                // size of each map in bytes
     
   uint64_t    oneParityMapIndex;      // used by the GPU
   uint64_t    dualParityMapM1Index;   // used by the GPU, not used for CisOneWithMultiSequenceHelper
   uint64_t    dualParityMapP1Index;   // used by the GPU, not used for CisOneWithMultiSequenceHelper
} v1_header_t;

class CisOneSequenceHelper : public AbstractSequenceHelper
{
public:
   CisOneSequenceHelper(App *theApp, uint64_t largestPrimeTested);

   ~CisOneSequenceHelper(void) {};

   void              CleanUp(void);
   
   void              LastChanceLogicBeforeSieving(void);
   
   void              ComputeSmallSieveLimit(void);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested) = 0;

   uint32_t          GetBaseMultiple(void) { return ii_BaseMultiple; };
   uint32_t          GetPowerResidueLcm(void) { return ii_PowerResidueLcm; };
   uint32_t          GetLimitBase(void) { return ii_LimitBase; };

   uint32_t          GetMaxBabySteps(void) { return ii_MaxBabySteps; };
   
   int16_t          *GetDivisorShifts(void) { return ip_DivisorShifts; };
   uint16_t         *GetPowerResidueIndices(void) { return ip_PowerResidueIndices; };
   
   legendre_t       *GetLegendre(void) { return ip_Legendre; };
   uint8_t          *GetLegendreTable(void) { return ip_LegendreTable; };
   
   void              BuildLegendreTables();

protected:
   void              BuildDivisorShifts(void);
   void              BuildPowerResidueIndices(void);
   virtual void      BuildCongruenceTables(void) = 0;
   bool              HasCongruentTerms(uint32_t ssIdx, uint32_t a, uint32_t b);
      
   virtual void      ComputeLegendreMemoryToAllocate(legendre_t *legendrePtr, uint64_t ssqfb) = 0;
   virtual void      AssignMemoryToLegendreTable(legendre_t *legendrePtr, uint64_t bytesUsed) = 0;
   virtual void      BuildLegendreTableForSequence(legendre_t *legendrePtr, uint64_t ssqfb, uint64_t stepsToDo, uint64_t stepsDone, time_t startTime) = 0;
   
   void              LoadLegendreTablesFromFile(legendre_t *legendrePtr);
   bool              ValidateLegendreFile(v1_header_t *headerPtr, legendre_t *legendrePtr);
   bool              ReadLegendreTableFromFile(FILE *fPtr, uint8_t *map, uint32_t mapSize, uint64_t offset);
   void              WriteLegendreTableToFile(legendre_t *legendrePtr);
   
   uint32_t          FindBestQ(uint32_t &expectedSubsequences);
   virtual double    RateQ(uint32_t Q, uint32_t s) = 0;

   // Set BASE_MULTIPLE to the smallest exponent Q for which sieving in
   // subsequence base b^Q will be considered. Must be a multiple of 2.
   uint32_t          ii_BaseMultiple;

   // For a prime p that satisfies p=1 (mod r), an "r-th power residue test"
   // checks whether a subsequence of k*b^n+c can possibly contain any terms of
   // the form x^r (mod p). If there are none then that subsequence can be
   // omitted from the BSGS step.
   //
   // To conduct r-th power residue tests for each r in a set R of prime powers,
   // set POWER_RESIDUE_LCM to lcm(R). POWER_RESIDUE_LCM must be a multiple of
   // BASE_MULTIPLE and must be less than 2^15.
   uint32_t          ii_PowerResidueLcm;

   // Allow sieving in base b^Q for Q chosen from the divisors of LIMIT_BASE.
   // Must be a multiple of POWER_RESIDUE_LCM.
   uint32_t          ii_LimitBase;
   
   uint32_t          ii_MaxBabySteps;

   int16_t          *ip_DivisorShifts;
   
   // This converts values r where POWER_RESIDUE_LCM % r == 0 to an index.  It is used to
   // reduce the memory requirements for ip_CongruentQIndices and ip_LadderIndices.   
   uint16_t         *ip_PowerResidueIndices;
   uint32_t          ii_UsedPowerResidueIndices;

   legendre_t       *ip_Legendre;
   uint8_t          *ip_LegendreTable;
   uint64_t          ii_LegendreBytes;

   inline uint64_t getNegCK(seq_t *seqPtr, uint64_t p)
   {
      uint64_t negCK;
      
      if (p < seqPtr->k)
         negCK = seqPtr->k % p;
      else 
         negCK = seqPtr->k;
     
      if (seqPtr->c > 0)
         negCK = p - negCK;
      
      return negCK;
   }
};

#endif

