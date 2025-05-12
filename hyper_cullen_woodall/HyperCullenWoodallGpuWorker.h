/* HyperCullenWoodallGpuWorker.h -- (C) Mark Rodenkirch

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _HyperCullenWoodallGpuWorker_H
#define _HyperCullenWoodallGpuWorker_H

#include "HyperCullenWoodallApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class HyperCullenWoodallGpuWorker : public Worker
{
public:
   HyperCullenWoodallGpuWorker(uint32_t myId, App *theApp);

   ~HyperCullenWoodallGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}  

   gpugroup_t       *ip_GpuGroup;
   uint64_t          il_NextTermsBuild;

   uint32_t          ii_GroupSize;
   uint32_t          ii_MaxGpuFactors;
   
   uint32_t         *ii_KernelTerms;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;

   HyperCullenWoodallApp          *ip_HyperCullenWoodallApp;

   GpuKernel        *ip_Kernel;
};

#endif

