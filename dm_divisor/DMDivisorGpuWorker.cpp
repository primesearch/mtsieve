/* DMDivisorWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "DMDivisorGpuWorker.h"
#include "dm_kernel.gpu.h"

DMDivisorGpuWorker::DMDivisorGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ip_DMDivisorApp = (DMDivisorApp *) theApp;
   
   ib_GpuWorker = true;
      
   il_MinK = ip_DMDivisorApp->GetMinK();
   il_MaxK = ip_DMDivisorApp->GetMaxK();
   
   ii_N = ip_DMDivisorApp->GetN();
      
   ii_MaxGpuFactors = ip_DMDivisorApp->GetMaxGpuFactors();
   
   snprintf(defines[defineCount++], 50, "#define THE_N %u\n", ii_N);
   snprintf(defines[defineCount++], 50, "#define K_MIN %" PRIu64"\n", il_MinK);
   snprintf(defines[defineCount++], 50, "#define K_MAX %" PRIu64"\n", il_MaxK);
   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("dm_kernel", dm_kernel, preKernelSources);

   ip_DMDivisorApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_DMDivisorApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 2*ii_MaxGpuFactors);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  DMDivisorGpuWorker::CleanUp(void)
{  
   delete ip_Kernel;
}

void  DMDivisorGpuWorker::TestMegaPrimeChunk(void)
{
   ii_FactorCount[0] = 0;
   
   ip_Kernel->Execute(ii_PrimesInList);

   for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
   {
      uint64_t k, prime;
      uint32_t idx = ii*2;
      
      k = (uint64_t) il_FactorList[idx+0];
      prime = il_FactorList[idx+1];
   
      ip_DMDivisorApp->ReportFactor(prime, k);
      
      if ((ii+1) == ii_MaxGpuFactors)
         break;
   }

   if (ii_FactorCount[0] >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0], ii_MaxGpuFactors);
      
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  DMDivisorGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("DMDivisorGpuWorker::TestMiniPrimeChunk not implemented");
}
