/* GenericGpuWorker.h -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GenericGpuWorker_H
#define _GenericGpuWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSubsequenceHelper.h"
#include "AbstractWorker.h"
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class GenericGpuWorker : public AbstractWorker
{
public:
   GenericGpuWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper);

   ~GenericGpuWorker() {};

   void              SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount);
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              InitializeWorker(void) {};
   
   uint32_t          ii_MaxGpuFactors;
   uint64_t         *il_FactorList;
   uint32_t          ii_FactorCount;
   
   uint64_t         *il_K;
   int64_t          *il_C;
   uint32_t         *ii_SeqIdx;
   uint32_t         *ii_Q;
   
   Kernel           *ip_SRKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_KASeqK;
   KernelArgument   *ip_KASeqC;
   KernelArgument   *ip_KASubSeqSeqIdx;
   KernelArgument   *ip_KASubSeqQ;
   KernelArgument   *ip_KAFactorCount;
   KernelArgument   *ip_KAFactorList;
};

#endif

