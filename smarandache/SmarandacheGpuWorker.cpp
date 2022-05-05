/* SmarandacheGpuWorker.cpp -- (C) Mark Rodenkirch, January 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "SmarandacheGpuWorker.h"
#include "sm_kernel.gpu.h"

SmarandacheGpuWorker::SmarandacheGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ib_GpuWorker = true;
   
   ip_SmarandacheApp = (SmarandacheApp *) theApp;

   ii_MaxGpuFactors = ip_SmarandacheApp->GetMaxGpuFactors();
   ii_MaxGpuSteps = ip_SmarandacheApp->GetMaxGpuSteps();
   
   terms_t *terms = ip_SmarandacheApp->GetTerms();
   
   ii_KernelCount = terms->termCount / ii_MaxGpuSteps;

   // In case it was rounded down
   if (ii_MaxGpuSteps * ii_KernelCount < terms->termCount)
      ii_KernelCount++;
   
   ii_AllTerms = (uint32_t **) xmalloc(ii_KernelCount * sizeof(uint32_t *));
   
   for (idx=0; idx<ii_KernelCount; idx++)
      ii_AllTerms[idx] = (uint32_t *) xmalloc((ii_MaxGpuSteps+1) * sizeof(uint32_t));
   
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = new Kernel(ip_SmarandacheApp->GetDevice(), "sm_kernel", sm_kernel, preKernelSources);

   ip_SmarandacheApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_SmarandacheApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   ii_KernelTerms = (uint32_t *) ip_Kernel->AddCpuArgument("factorCount", sizeof(uint32_t), ii_MaxGpuSteps+1);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 2*ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(0);
   
   ib_NeedToRebuildTerms = true;
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  SmarandacheGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
   
   for (uint32_t idx=0; idx<ii_KernelCount; idx++)
      xfree(ii_AllTerms[idx]);
}

void  SmarandacheGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t kIdx;
   uint32_t n, ii, idx;
   uint64_t prime;
   time_t   reportTime;

   if (ib_NeedToRebuildTerms)
      BuildTerms(ip_SmarandacheApp->GetTerms());

   for (kIdx=0; kIdx<ii_KernelCount; kIdx++)
   {
      // Check if we need to execute this kernel
      if (ii_AllTerms[kIdx][0] == 0)
         break;

      memcpy(ii_KernelTerms, ii_AllTerms[kIdx], ii_MaxGpuSteps * sizeof(uint32_t));
      
      reportTime = time(NULL) + 60;

      ii_FactorCount[0] = 0;
  
      ip_Kernel->Execute(ii_WorkSize);

      for (ii=0; ii<ii_FactorCount[0]; ii++)
      {  
         idx = ii*2;
         
         n = (uint32_t) il_FactorList[idx+0];
         prime = il_FactorList[idx+1];
      
         ip_SmarandacheApp->ReportFactor(prime, n);
         
         ib_NeedToRebuildTerms = true;
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0], ii_MaxGpuFactors);

      if (kIdx < ii_KernelCount && ip_SmarandacheApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_SmarandacheApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, kIdx, ii_KernelCount);
         reportTime = time(NULL) + 60;
      }
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  SmarandacheGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheGpuWorker::TestMiniPrimeChunk not implemented");
}

void  SmarandacheGpuWorker::BuildTerms(terms_t *terms)
{
   uint32_t kernelIdx, kernelTermsIdx, globalTermsIdx;
   uint32_t maxGpuSteps = ip_SmarandacheApp->GetMaxGpuSteps();
   
   globalTermsIdx = 0;
   for (kernelIdx=0; kernelIdx<ii_KernelCount; kernelIdx++)
   {
      memset(ii_AllTerms[kernelIdx], 0x00, (maxGpuSteps + 1) * sizeof(uint32_t));
      
      if (globalTermsIdx >= terms->termCount)
         continue;

      for (kernelTermsIdx=0; kernelTermsIdx<maxGpuSteps; kernelTermsIdx++)
      {
         ii_AllTerms[kernelIdx][kernelTermsIdx] = terms->termList[globalTermsIdx];
         globalTermsIdx++;
         
         if (globalTermsIdx >= terms->termCount)
            break;
      };
   }
   
   ib_NeedToRebuildTerms = false;
}
