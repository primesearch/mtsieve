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
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

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
   uint64_t         *ip_FactorList;
   uint32_t          ii_FactorCount;
   
   uint32_t         *ip_BabySteps;
   uint32_t         *ip_GiantSteps;
      
   uint32_t         *ip_CongruentQIndices;
   uint32_t         *ip_LadderIndices;
   
   uint16_t         *ip_AllQs;
   uint16_t         *ip_AllLadders;
   
   uint16_t         *ip_SeqQs;
   uint16_t         *ip_SeqLadders;
   
   legendre_t       *ip_Legendre;
   uint8_t          *ip_LegendreTable;
   
   uint32_t          ii_Dim1;
   uint32_t          ii_Dim2;
   uint32_t          ii_Dim3;
   
   Kernel           *ip_SRKernel;
   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_KADualParityMapM1;
   KernelArgument   *ip_KADualParityMapP1;
   KernelArgument   *ip_KAOneParityMap;
   KernelArgument   *ip_KASubSeqBS;
   KernelArgument   *ip_KASubSeqGS;
   KernelArgument   *ip_KADivisorShifts;
   KernelArgument   *ip_KAPrlIndices;
   KernelArgument   *ip_KAQIndices;
   KernelArgument   *ip_KAQs;
   KernelArgument   *ip_KALadderIndices;
   KernelArgument   *ip_KALadders;
   KernelArgument   *ip_KAFactorCount;
   KernelArgument   *ip_KAFactorList;
};

#endif

