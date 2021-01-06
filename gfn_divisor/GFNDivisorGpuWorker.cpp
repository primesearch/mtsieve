/* GFNDivisorWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "GFNDivisorGpuWorker.h"
#include "gfn_kernel.h"

GFNDivisorGpuWorker::GFNDivisorGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{    
   const char *source[10];
   char        define1[50];
   char        define2[50];
   char        define3[50];
   char        define4[50];
   char        define5[50];
   
   ip_GFNDivisorApp = (GFNDivisorApp *) theApp;
   
   ib_GpuWorker = true;
      
   il_MinK = ip_GFNDivisorApp->GetMinK();
   il_MaxK = ip_GFNDivisorApp->GetMaxK();
   
   ii_MinN = ip_GFNDivisorApp->GetMinN();
   ii_MaxN = ip_GFNDivisorApp->GetMaxN();
      
   ii_MaxGpuFactors = ip_GFNDivisorApp->GetMaxGpuFactors();
   
   sprintf(define1, "#define N_MIN %u\n", ii_MinN);
   sprintf(define2, "#define N_MAX %u\n", ii_MaxN);
   sprintf(define3, "#define K_MIN %" PRIu64"\n", il_MinK);
   sprintf(define4, "#define K_MAX %" PRIu64"\n", il_MaxK);
   sprintf(define5, "#define MAX_FACTORS %d\n", ii_MaxGpuFactors);
   
   source[0] = define1;
   source[1] = define2;
   source[2] = define3;
   source[3] = define4;
   source[4] = define5;
   source[5] = gfn_kernel;
   source[6] = 0;

   ip_GFNDivisorKernel = new Kernel(ip_GFNDivisorApp->GetDevice(), "gfn_kernel", source);

   AllocatePrimeList(ip_GFNDivisorKernel->GetWorkGroupSize());

   il_FactorList  =  (uint64_t *) xmalloc(4*ii_MaxGpuFactors*sizeof(uint64_t));

   ip_KAPrime        = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KAFactorCount  = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList   = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);

   ip_GFNDivisorKernel->AddArgument(ip_KAPrime);
   ip_GFNDivisorKernel->AddArgument(ip_KAFactorCount);
   ip_GFNDivisorKernel->AddArgument(ip_KAFactorList);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  GFNDivisorGpuWorker::CleanUp(void)
{  
   delete ip_KAPrime;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;
   delete ip_GFNDivisorKernel;

   xfree(il_FactorList);
}

void  GFNDivisorGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, n;
   uint64_t k, prime;
   
   ii_FactorCount = 0;
   
   ip_GFNDivisorKernel->Execute(ii_WorkSize);

   for (uint32_t ii=0; ii<ii_FactorCount; ii++)
   {  
      idx = ii*4;
      
      k = (uint64_t) il_FactorList[idx+0];
      n = (uint32_t) il_FactorList[idx+1];
      prime = il_FactorList[idx+2];
   
      ip_GFNDivisorApp->ReportFactor(prime, k, n, true);
      
      if (ii_FactorCount >= ii_MaxGpuFactors)
         break;
   }

   if (ii_FactorCount >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount, ii_MaxGpuFactors);
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  GFNDivisorGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GFNDivisorGpuWorker::TestMiniPrimeChunk not implemented");
}
