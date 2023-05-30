/* XYYXSparseGpuWorker.h -- (C) Mark Rodenkirch, June 2014

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _XYYXSparseGpuWorker_H
#define _XYYXSparseGpuWorker_H

#include "XYYXApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class XYYXSparseGpuWorker : public Worker
{
public:
   XYYXSparseGpuWorker(uint32_t myId, App *theApp);

   ~XYYXSparseGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}  
   gputerm_t        *ip_Terms;
   uint64_t          il_NextTermsBuild;
   
   uint32_t          ii_Groups;
   uint32_t          ii_GroupSize;
   uint32_t          ii_MaxGpuFactors;
   
   uint32_t         *ii_KernelTerms;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;

   XYYXApp          *ip_XYYXApp;

   GpuKernel        *ip_Kernel;
};

#endif

