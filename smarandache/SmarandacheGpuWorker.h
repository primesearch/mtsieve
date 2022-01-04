/* SmarandacheGpuWorker.h -- (C) Mark Rodenkirch, January 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SmarandacheGpuWorker_H
#define _SmarandacheGpuWorker_H

#include "SmarandacheApp.h"
#include "../core/Worker.h"
#include "../opencl/Kernel.h"
#include "../opencl/KernelArgument.h"

using namespace std;

class SmarandacheGpuWorker : public Worker
{
public:
   SmarandacheGpuWorker(uint32_t myId, App *theApp);

   ~SmarandacheGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   SmarandacheApp *ip_SmarandacheApp;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_Smarandache;
   uint32_t          ii_MaxGpuFactors;
   
   uint64_t         *il_RemainderList;
   uint32_t          ii_Params[5];
   uint32_t          ii_FactorCount;
   int64_t          *il_FactorList;

   Kernel           *ip_FactorialKernel;

   KernelArgument   *ip_KAPrime;
   KernelArgument   *ip_KARemainder;
   KernelArgument   *ip_KAParams;
   KernelArgument   *ip_KATemp;
   KernelArgument   *ip_KAFactorCount;
   KernelArgument   *ip_KAFactorList;
};

#endif

