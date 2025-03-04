/* SmarandacheWellinGpuWorker.h -- (C) Mark Rodenkirch, March 2025

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SmarandacheWellinGpuWorker_H
#define _SmarandacheWellinGpuWorker_H

#include "SmarandacheWellinApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class SmarandacheWellinGpuWorker : public Worker
{
public:
   SmarandacheWellinGpuWorker(uint32_t myId, App *theApp);

   ~SmarandacheWellinGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}
   
   SmarandacheWellinApp   *ip_SmarandacheWellinApp;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   uint64_t          il_NextTermsBuild;
      
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_MaxGpuSteps;
   
   uint32_t         *ii_KernelTerms;
   uint16_t         *ii_KernelTermGaps;
   uint8_t          *ii_RemainingTerms;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;

   GpuKernel        *ip_Kernel;
};

#endif

