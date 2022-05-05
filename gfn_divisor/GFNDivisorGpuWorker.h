/* GFNDivisorWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GFNDivisorGpuWorker_H
#define _GFNDivisorGpuWorker_H

#include "GFNDivisorApp.h"
#include "../core/Worker.h"

#ifdef USE_OPENCL
#include "../gpu_opencl/Kernel.h"
#else
#include "../gpu_metal/Kernel.h"
#endif


using namespace std;

class GFNDivisorGpuWorker : public Worker
{
public:
   GFNDivisorGpuWorker(uint32_t myId, App *theApp);

   ~GFNDivisorGpuWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   GFNDivisorApp    *ip_GFNDivisorApp;
      
   Kernel           *ip_Kernel;
   
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_MaxIterations;
   
   uint32_t         *ii_Parameters;
   uint64_t         *il_Remainders;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;
};

#endif

