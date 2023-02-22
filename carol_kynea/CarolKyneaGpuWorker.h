/* CarolKyneaGpuWorker.h -- (C) Mark Rodenkirch, November 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CarolKyneaGpuWorker_H
#define _CarolKyneaGpuWorker_H

#include "CarolKyneaApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class CarolKyneaGpuWorker : public Worker
{
public:
   CarolKyneaGpuWorker(uint32_t myId, App *theApp);

   ~CarolKyneaGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}
   uint32_t          ii_MaxGpuFactors;

   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;

   // Variables used by the kernel
   uint32_t         *ii_FactorCount;
   int64_t          *il_FactorList;

   CarolKyneaApp *ip_CarolKyneaApp;
   
   GpuKernel        *ip_Kernel;
};

#endif

