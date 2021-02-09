/* PrimorialGpuWorker.h -- (C) Mark Rodenkirch, September 2016

   This is the header for the class sets up the call to the fsieve GPU function
   and parses the output from the GPU to determine if we have a factor.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimorialGpuWorker_H
#define _PrimorialGpuWorker_H

#include "PrimorialApp.h"
#include "../core/Worker.h"
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class PrimorialGpuWorker : public Worker
{
public:
   PrimorialGpuWorker(uint32_t myId, App *theApp);

   ~PrimorialGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   PrimorialApp *ip_PrimorialApp;
   
   uint32_t          ii_MinPrimorial;
   uint32_t          ii_MaxPrimorial;
   
   uint32_t         *ip_PrimorialPrimes;
   uint32_t          ii_NumberOfPrimorialPrimes;
   
   uint16_t         *ip_PrimorialPrimeGaps;
   uint16_t          ii_BiggestGap;
   
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
   
   uint64_t         *ip_RemainderList;
   uint64_t          ii_Params[5];
   uint32_t          ii_MaxIterations;
   uint32_t          ii_FactorCount;
   int64_t          *ip_FactorList;

   Kernel           *ip_PrimorialKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_KAPrimorialGaps;
   KernelArgument   *ip_KARemainder;
   KernelArgument   *ip_KAParams;
   KernelArgument   *ip_KAFactorCount;
   KernelArgument   *ip_KAFactorList;
};

#endif

