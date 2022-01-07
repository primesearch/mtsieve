/* SmarandacheGpuWorker.cpp -- (C) Mark Rodenkirch, January 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "SmarandacheGpuWorker.h"
#include "sm_kernel.h"
#include <time.h>

SmarandacheGpuWorker::SmarandacheGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        define1[40];
   const char *source[10];
   uint32_t    maxGpuSteps, idx;

   ib_GpuWorker = true;
   
   ip_SmarandacheApp = (SmarandacheApp *) theApp;

   ii_MaxGpuFactors = ip_SmarandacheApp->GetMaxGpuFactors();
   maxGpuSteps = ip_SmarandacheApp->GetMaxGpuSteps();
   
   terms_t *terms = ip_SmarandacheApp->GetTerms();
   
   ii_KernelCount = terms->termCount / maxGpuSteps;

   // In case it was rounded down
   if (maxGpuSteps * ii_KernelCount < terms->termCount)
      ii_KernelCount++;
   
   ii_Terms = (uint32_t **) xmalloc(ii_KernelCount * sizeof(uint32_t *));
   ip_SmarandacheKernel = (Kernel **) xmalloc(ii_KernelCount * sizeof(Kernel *));
   ip_KATerms = (KernelArgument **) xmalloc(ii_KernelCount * sizeof(KernelArgument *));
   
   for (idx=0; idx<ii_KernelCount; idx++)
      ii_Terms[idx] = (uint32_t *) xmalloc((maxGpuSteps+1) * sizeof(uint32_t));

   BuildTerms(terms);
   
   sprintf(define1, "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   
   source[0] = define1;
   source[1] = sm_kernel;
   source[2] = 0;

   ip_SmarandacheKernel[0] = new Kernel(ip_SmarandacheApp->GetDevice(), "sm_kernel", source);
      
   AllocatePrimeList(ip_SmarandacheKernel[0]->GetWorkGroupSize());
   
   delete ip_SmarandacheKernel[0];
   
   il_FactorList    = (uint64_t *)  xmalloc(2*ii_MaxGpuFactors*sizeof(uint64_t));

   ip_KAPrime        = new KernelArgument(ip_SmarandacheApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KAFactorCount  = new KernelArgument(ip_SmarandacheApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList   = new KernelArgument(ip_SmarandacheApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 2*ii_MaxGpuFactors);

   for (idx=0; idx<ii_KernelCount; idx++)
   {
      ip_KATerms[idx]  = new KernelArgument(ip_SmarandacheApp->GetDevice(), "terms", KA_HOST_TO_GPU, ii_Terms[idx], maxGpuSteps+1);
      
      ip_SmarandacheKernel[idx] = new Kernel(ip_SmarandacheApp->GetDevice(), "sm_kernel", source);

      ip_SmarandacheKernel[idx]->AddArgument(ip_KAPrime);
      ip_SmarandacheKernel[idx]->AddArgument(ip_KATerms[idx]);
      ip_SmarandacheKernel[idx]->AddArgument(ip_KAFactorCount);
      ip_SmarandacheKernel[idx]->AddArgument(ip_KAFactorList);
   }
   
   ip_SmarandacheKernel[0]->PrintStatistics(0);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  SmarandacheGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;
   
   for (uint32_t idx=0; idx<ii_KernelCount; idx++)
   {
      delete ip_SmarandacheKernel[idx];
      delete ip_KATerms[idx];
      xfree(ii_Terms[idx]);
   }
   
   xfree(il_PrimeList);
   xfree(il_FactorList);
   
   xfree(ip_SmarandacheKernel);
}

void  SmarandacheGpuWorker::BuildTerms(terms_t *terms)
{
   uint32_t kernelIdx, kernelTermsIdx, globalTermsIdx;
   uint32_t maxGpuSteps = ip_SmarandacheApp->GetMaxGpuSteps();
   
   globalTermsIdx = 0;
   for (kernelIdx=0; kernelIdx<ii_KernelCount; kernelIdx++)
   {
      memset(ii_Terms[kernelIdx], 0x00, (maxGpuSteps + 1) * sizeof(uint32_t));
      
      if (globalTermsIdx >= terms->termCount)
         continue;

      for (kernelTermsIdx=0; kernelTermsIdx<maxGpuSteps; kernelTermsIdx++)
      {
         ii_Terms[kernelIdx][kernelTermsIdx] = terms->termList[globalTermsIdx];
         globalTermsIdx++;
         
         if (globalTermsIdx >= terms->termCount)
            break;
      };
   }
   
   ib_NeedToRebuild = false;
}

void  SmarandacheGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t kernelIdx;
   uint32_t n, ii, idx;
   uint64_t prime;
   time_t   reportTime;

   if (ib_NeedToRebuild)
      BuildTerms(ip_SmarandacheApp->GetTerms());

   for (kernelIdx=0; kernelIdx<ii_KernelCount; kernelIdx++)
   {
      // Check if we need to execute this kernel
      if (ii_Terms[kernelIdx][0] == 0)
         break;

      reportTime = time(NULL) + 60;

      ii_FactorCount = 0;
  
      ip_SmarandacheKernel[kernelIdx]->Execute(ii_WorkSize);

      for (ii=0; ii<ii_FactorCount; ii++)
      {  
         idx = ii*2;
         
         n = (uint32_t) il_FactorList[idx+0];
         prime = il_FactorList[idx+1];
      
         ip_SmarandacheApp->ReportFactor(prime, n);
         
         ib_NeedToRebuild = true;
         
         if (ii >= ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors.  Use -M to increase max factors", ii_FactorCount);

      if (kernelIdx < ii_KernelCount && ip_SmarandacheApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_SmarandacheApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, kernelIdx, ii_KernelCount);
         reportTime = time(NULL) + 60;
      }
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  SmarandacheGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheGpuWorker::TestMiniPrimeChunk not implemented");
}
