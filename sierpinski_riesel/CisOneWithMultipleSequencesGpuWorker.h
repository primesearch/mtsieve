/* CisOneWithMultipleSequencesGpuWorker.h -- (C) Mark Rodenkirch, July 2023

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CisOneWithMultipleSequencesGpuWorker_H
#define _CisOneWithMultipleSequencesGpuWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSequenceHelper.h"
#include "AbstractWorker.h"
#include "CisOneWithMultipleSequencesHelper.h"

#include "../core/GpuKernel.h"

class CisOneWithMultipleSequencesGpuWorker : public AbstractWorker
{
public:
   CisOneWithMultipleSequencesGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper);

   ~CisOneWithMultipleSequencesGpuWorker() {};

   void              Prepare(uint64_t largestPrimeTested, uint32_t bestQ);
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {};
   void              CreateKernels(uint32_t sequencesPerKernel);
   void              PopulateKernelArguments(void);
   void              CreateKernel(void);

   uint32_t          ii_KernelCount;
   uint32_t          ii_SequencesPerKernel;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_MaxSequences;
   uint32_t          ii_MaxSubsequences;

   CisOneWithMultipleSequencesHelper *ip_CisOneHelper;   
   
   GpuKernel        *ip_Kernel;

   uint64_t         *il_Primes;
   
   uint64_t         *il_Ks;
   int64_t          *il_KCCores;
   uint32_t         *ii_SeqData;
   uint32_t         *ii_SubseqData;
                    
   int16_t          *ii_DivisorShifts;
   uint16_t         *ii_PowerResidueIndices;
   uint32_t         *ii_CongruentSubseqIndices;
   uint32_t         *ii_CongruentSubseqs;
                    
   uint16_t         *ii_Ladders;

   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;
};

#endif

