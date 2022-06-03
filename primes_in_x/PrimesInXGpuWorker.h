/* PixSieve.h -- (C) Mark Rodenkirch, September 2016

   This is the header for the class sets up the call to the fsieve GPU function
   and parses the output from the GPU to determine if we have a factor.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _PrimesInXGpuWorker_H
#define _PrimesInXGpuWorker_H

#include "PrimesInXApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class PrimesInXGpuWorker : public Worker
{
public:
   PrimesInXGpuWorker(uint32_t myId, App *theApp);

   ~PrimesInXGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   uint32_t         *BuildDigitList(void);
   uint32_t         *BuildDigitList(uint32_t multiplier, uint32_t power, uint32_t *digitList);

   uint32_t          ii_MaxSteps;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_GroupSize;
   
   uint32_t         *ii_e1DigitList;
   uint32_t         *ii_e3DigitList;
   uint32_t         *ii_e6DigitList;
   uint32_t         *ii_e9DigitList;
   
   uint64_t         *il_Residuals;
   uint32_t         *ii_DigitList;
   uint32_t         *ii_FactorCount;
   uint64_t         *il_FactorList;

   PrimesInXApp     *ip_PrimesInXApp;

   GpuKernel        *ip_Kernel;
};

#endif

