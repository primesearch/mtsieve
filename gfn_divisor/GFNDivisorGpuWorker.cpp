/* GFNDivisorWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "GFNDivisorGpuWorker.h"
#include "gfn_kernel.gpu.h"

GFNDivisorGpuWorker::GFNDivisorGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ip_GFNDivisorApp = (GFNDivisorApp *) theApp;
   
   ib_GpuWorker = true;
      
   il_MinK = ip_GFNDivisorApp->GetMinK();
   il_MaxK = ip_GFNDivisorApp->GetMaxK();
   
   ii_MinN = ip_GFNDivisorApp->GetMinN();
   ii_MaxN = ip_GFNDivisorApp->GetMaxN();
      
   ii_MaxGpuSteps = ip_GFNDivisorApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_GFNDivisorApp->GetMaxGpuFactors();
   
   ii_MaxIterations = 1 + ((1 + ii_MaxN - ii_MinN) / ii_MaxGpuSteps);

   sprintf(defines[defineCount++], "#define N_MIN %u\n", ii_MinN);
   sprintf(defines[defineCount++], "#define N_MAX %u\n", ii_MaxN);
   sprintf(defines[defineCount++], "#define K_MIN %" PRIu64"\n", il_MinK);
   sprintf(defines[defineCount++], "#define K_MAX %" PRIu64"\n", il_MaxK);
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   
   if (ii_MaxIterations > 1)
   {
      sprintf(defines[defineCount++], "#define D_MULTI_PASS\n");
      sprintf(defines[defineCount++], "#define D_MAX_STEPS %d\n", ii_MaxGpuSteps);
   }
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("gfn_kernel", gfn_kernel, preKernelSources);

   ip_GFNDivisorApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_GFNDivisorApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   
   if (ii_MaxIterations > 1)
   {
      ii_Parameters = (uint32_t *) ip_Kernel->AddCpuArgument("parameters", sizeof(uint32_t), 5);
      il_Remainders = (uint64_t *) ip_Kernel->AddSharedArgument("remainders", sizeof(uint64_t), ii_PrimesInList);
   }
      
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  GFNDivisorGpuWorker::CleanUp(void)
{  
   delete ip_Kernel;
}

void  GFNDivisorGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, n, iteration;
   uint64_t k, prime;

   for (iteration=0; iteration<ii_MaxIterations; iteration++)
   {
      if (ii_MaxIterations > 1)
         ii_Parameters[0] = ii_MinN + (iteration * ii_MaxGpuSteps);
      
      ii_FactorCount[0] = 0;
      
      ip_Kernel->Execute(ii_PrimesInList);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {
         idx = ii*4;
         
         k = (uint64_t) il_FactorList[idx+0];
         n = (uint32_t) il_FactorList[idx+1];
         prime = il_FactorList[idx+2];
      
         ip_GFNDivisorApp->ReportFactor(prime, k, n, true);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0], ii_MaxGpuFactors);
   }
      
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  GFNDivisorGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GFNDivisorGpuWorker::TestMiniPrimeChunk not implemented");
}
