/* CisOneWithOneSequenceGpuWorker.h -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CisOneWithOneSequenceGpuWorker_H
#define _CisOneWithOneSequenceGpuWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSequenceHelper.h"
#include "CisOneWithOneSequenceHelper.h"
#include "AbstractWorker.h"

#include "../core/GpuKernel.h"

using namespace std;

class CisOneWithOneSequenceGpuWorker : public AbstractWorker
{
public:
   CisOneWithOneSequenceGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper);

   ~CisOneWithOneSequenceGpuWorker() {};

   void              Prepare(uint64_t largestPrimeTested, uint32_t bestQ);
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   CisOneWithOneSequenceHelper *ip_CisOneHelper;
   
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_ChunksPerGpuWorker;
   uint32_t          ii_KernelWorkSize;
   
   GpuKernel        *ip_Kernel;
   
   uint64_t         *il_Primes;
   uint8_t          *ii_DualParityMapM1;
   uint8_t          *ii_DualParityMapP1;
   uint8_t          *ii_SingleParityMap;
   
   uint32_t         *ii_BabySteps;
   uint32_t         *ii_GiantSteps;

   int16_t          *ii_DivisorShifts;
   uint16_t         *ii_PowerResidueIndices;

   uint32_t         *ii_QIndices;
   uint16_t         *ii_Qs;
   
   uint32_t         *ii_LadderIndices;
   uint16_t         *ii_Ladders;
   
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;
   
   
};

#endif

