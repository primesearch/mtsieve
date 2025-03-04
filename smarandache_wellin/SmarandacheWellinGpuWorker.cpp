/* SmarandacheWellinGpuWorker.cpp -- (C) Mark Rodenkirch, March 2025

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "SmarandacheWellinGpuWorker.h"
#include "smw_kernel.gpu.h"

SmarandacheWellinGpuWorker::SmarandacheWellinGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   uint32_t    numberOfPrimes;
   uint16_t    biggestGap;
   
   ib_GpuWorker = true;
   
   ip_SmarandacheWellinApp = (SmarandacheWellinApp *) theApp;

   ii_MaxGpuFactors = ip_SmarandacheWellinApp->GetMaxGpuFactors();
   ii_MaxGpuSteps = ip_SmarandacheWellinApp->GetMaxGpuSteps();
   
   ip_SmarandacheWellinApp = (SmarandacheWellinApp *) theApp;
   
   ii_MinN = ip_SmarandacheWellinApp->GetMinN();
   ii_MaxN = ip_SmarandacheWellinApp->GetMaxN();

   uint32_t *primes = ip_SmarandacheWellinApp->GetPrimes(numberOfPrimes);
   uint16_t *primeGaps = ip_SmarandacheWellinApp->GetPrimeGaps(biggestGap);

   snprintf(defines[defineCount++], 50, "#define D_MIN_N %d\n", ii_MinN);
   snprintf(defines[defineCount++], 50, "#define D_MAX_N %d\n", ii_MaxN);
   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   snprintf(defines[defineCount++], 50, "#define D_BIGGEST_GAP %d\n", biggestGap);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("smw_kernel", smw_kernel, preKernelSources);

   ip_SmarandacheWellinApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_SmarandacheWellinApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   ii_KernelTerms = (uint32_t *) ip_Kernel->AddCpuArgument("terms", sizeof(uint32_t), numberOfPrimes+1);
   ii_KernelTermGaps = (uint16_t *) ip_Kernel->AddCpuArgument("termGaps", sizeof(uint16_t), numberOfPrimes+1);
   ii_RemainingTerms = (uint8_t *) ip_Kernel->AddCpuArgument("remainingTerms", sizeof(uint8_t), numberOfPrimes+1);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 2*ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(0);

   for (uint32_t i=0; i<numberOfPrimes; i++)
   {
      ii_KernelTerms[i] = primes[i];
      ii_KernelTermGaps[i] = primeGaps[i];
   }
   
   il_NextTermsBuild = 0;
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  SmarandacheWellinGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  SmarandacheWellinGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t n, ii, idx;
   uint64_t prime;
   
   if (il_NextTermsBuild < il_PrimeList[0])
   {
      ip_SmarandacheWellinApp->FillTerms(ii_RemainingTerms);
      
      il_NextTermsBuild = il_PrimeList[0] * 2;
   }

   ii_FactorCount[0] = 0;

   ip_Kernel->Execute(ii_PrimesInList);

   for (ii=0; ii<ii_FactorCount[0]; ii++)
   {  
      idx = ii*2;
      
      n = (uint32_t) il_FactorList[idx+0];
      prime = il_FactorList[idx+1];
      
      ip_SmarandacheWellinApp->ReportFactor(prime, n);
      
      if ((ii+1) == ii_MaxGpuFactors)
         break;
   }

   if (ii_FactorCount[0] >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0], ii_MaxGpuFactors);

   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  SmarandacheWellinGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheWellinGpuWorker::TestMiniPrimeChunk not implemented");
}
