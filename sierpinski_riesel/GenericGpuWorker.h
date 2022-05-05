/* GenericGpuWorker.h -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GenericGpuWorker_H
#define _GenericGpuWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSequenceHelper.h"
#include "AbstractWorker.h"

#ifdef USE_OPENCL
#include "../gpu_opencl/Kernel.h"
#else
#include "../gpu_metal/Kernel.h"
#endif

using namespace std;

class GenericGpuWorker : public AbstractWorker
{
public:
   GenericGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper);

   ~GenericGpuWorker() {};

   void              Prepare(uint64_t largestPrimeTested, uint32_t bestQ);
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   Kernel           *CreateKernel(uint32_t kernelIdx, uint32_t sequences, uint32_t subsequences);
   
   bool              ib_CanUseCIsOneLogic;
   uint64_t          il_MaxK;
   
   uint32_t          ii_KernelCount;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_SequencesPerKernel;
   
   uint32_t         *ii_SubseqIdx;
   
   uint64_t        **il_Primes;
   uint64_t        **il_K;
   int64_t         **il_C;
   uint32_t        **ii_SeqIdx;
   uint32_t        **ii_Q;
   uint32_t        **ii_FactorCount;
   uint64_t        **il_FactorList;
   
   Kernel          **ip_Kernel;
};

#endif

