/* LifchitzGpuWorker.cpp -- (C) Mark Rodenkirch, April 2024

   This class sets up the call to the LifchitzGpuWorker GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "LifchitzGpuWorker.h"
#include "lifchitz_kernel.gpu.h"

#define HASH_MAX_DENSITY    0.65
#define HASH_MINIMUM_ELTS   8
#define HASH_MINIMUM_SHIFT 11

LifchitzGpuWorker::LifchitzGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   char        defines[20][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ib_GpuWorker = true;
   
   ip_LifchitzApp = (LifchitzApp *) theApp;

   ii_MinX = ip_LifchitzApp->GetMinX();
   ii_MaxX = ip_LifchitzApp->GetMaxX();
   ii_MinY = ip_LifchitzApp->GetMinY();
   ii_MaxY = ip_LifchitzApp->GetMaxY();

   ii_XPerChunk = ip_LifchitzApp->GetXPerChunk();
   
   if (ii_XPerChunk > ii_MaxX - ii_MinX)
      ii_XPerChunk = 1 + ii_MaxX - ii_MinX;
   
   ii_MaxGpuFactors = ip_LifchitzApp->GetMaxGpuFactors();

   // This allows us to have
   uint32_t hashElements = ii_XPerChunk;
   
   if (hashElements < HASH_MINIMUM_ELTS)
      hashElements = HASH_MINIMUM_ELTS;
   
   uint32_t hashSize = (1 << HASH_MINIMUM_SHIFT);
   
   while (hashSize < hashElements/HASH_MAX_DENSITY)
      hashSize <<= 1;

   snprintf(defines[defineCount++], 50, "#define D_MIN_X %u", ii_MinX);
   snprintf(defines[defineCount++], 50, "#define D_MAX_X %u", ii_MaxX);
   snprintf(defines[defineCount++], 50, "#define D_MIN_Y %u", ii_MinY);
   snprintf(defines[defineCount++], 50, "#define D_MAX_Y %u", ii_MaxY);
   snprintf(defines[defineCount++], 50, "#define X_PER_CHUNK %u", ii_XPerChunk);
   snprintf(defines[defineCount++], 50, "#define HASH_ELEMENTS %u", hashElements);
   snprintf(defines[defineCount++], 50, "#define HASH_SIZE %u", hashSize);
   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("lifchitz_kernel", lifchitz_kernel, preKernelSources);

   ip_LifchitzApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_LifchitzApp->GetGpuPrimesPerWorker();
   
   il_PrimeList    = (uint64_t *)  ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t),  ii_PrimesInList);
   ii_FactorCount  = (uint32_t *)  ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 10);
   ip_FactorList   = (factor_t *)  ip_Kernel->AddGpuArgument("factorList", sizeof(factor_t), ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(0);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  LifchitzGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  LifchitzGpuWorker::TestMegaPrimeChunk(void)
{
   ii_FactorCount[0] = 0;
   
   ip_Kernel->Execute(ii_PrimesInList);
   
   for (uint32_t fIdx=0; fIdx<ii_FactorCount[0]; fIdx++)
   {
      if (fIdx == ii_MaxGpuFactors)
         break;
      
      uint32_t x = (uint32_t) ip_FactorList[fIdx].x;
      uint32_t y = (uint32_t) ip_FactorList[fIdx].y;
      int32_t  sign = (ip_FactorList[fIdx].sign == 1 ? +1 : -1);
      uint64_t theFactor = ip_FactorList[fIdx].factor;
            
      ip_LifchitzApp->ReportFactor(theFactor, x, y, sign);
   }

   if (ii_FactorCount[0] >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
 
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  LifchitzGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("LifchitzGpuWorker::TestMiniPrimeChunk not implemented");
}
