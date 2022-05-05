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

#ifdef USE_OPENCL
#include "../gpu_opencl/Kernel.h"
#else
#include "../gpu_metal/Kernel.h"
#endif

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
   void              BuildTerms(terms_t *terms);
   
   SmarandacheApp   *ip_SmarandacheApp;

   bool              ib_NeedToRebuildTerms;
   uint32_t          ii_KernelCount;
   uint32_t        **ii_AllTerms;
   
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_MaxGpuSteps;
   
   uint32_t         *ii_KernelTerms;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;

   Kernel           *ip_Kernel;
};

#endif

