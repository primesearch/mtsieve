/* DMDivisorWorker.h -- (C) Mark Rodenkirch, October 2025

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _DMDivisorGpuWorker_H
#define _DMDivisorGpuWorker_H

#include "DMDivisorApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class DMDivisorGpuWorker : public Worker
{
public:
   DMDivisorGpuWorker(uint32_t myId, App *theApp);

   ~DMDivisorGpuWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList) {}

private:
   DMDivisorApp     *ip_DMDivisorApp;
      
   GpuKernel        *ip_Kernel;
   
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_N;
   uint32_t          ii_MaxGpuFactors;
   
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;
};

#endif

