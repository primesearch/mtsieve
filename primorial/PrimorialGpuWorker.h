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

#include "../core/GpuKernel.h"

class PrimorialGpuWorker : public Worker
{
public:
   PrimorialGpuWorker(uint32_t myId, App *theApp);

   ~PrimorialGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              BuildPrimeGapGroups(void);
   
   PrimorialApp     *ip_PrimorialApp;
   
   uint32_t          ii_MinPrimorial;
   uint32_t          ii_MaxPrimorial;
      
   uint32_t          ii_IdxOfMinPrimorial;
   uint32_t          ii_GroupCount;
   uint16_t          ii_BiggestGap;
   uint16_t        **ip_GapGroups;
   
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MaxGpuFactors;
   
   uint64_t         *il_Residues;
   uint16_t         *ii_PrimeGaps;
   uint32_t         *ii_Parameters;
   uint32_t         *ii_FactorCount;
   int64_t          *il_FactorList;

   GpuKernel        *ip_Kernel;
};

#endif

