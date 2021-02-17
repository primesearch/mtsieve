/* GFNDivisorWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
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
   char        define6[50];
   char        define7[50];
   
   ip_GFNDivisorApp = (GFNDivisorApp *) theApp;
   
   ib_GpuWorker = true;
      
   il_MinK = ip_GFNDivisorApp->GetMinK();
   il_MaxK = ip_GFNDivisorApp->GetMaxK();
   
   ii_MinN = ip_GFNDivisorApp->GetMinN();
   ii_MaxN = ip_GFNDivisorApp->GetMaxN();
      
   ii_MaxGpuSteps = ip_GFNDivisorApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_GFNDivisorApp->GetMaxGpuFactors();
   
   ii_MaxIterations = 1 + ((1 + ii_MaxN - ii_MinN) / ii_MaxGpuSteps);
   
   sprintf(define1, "#define N_MIN %u\n", ii_MinN);
   sprintf(define2, "#define N_MAX %u\n", ii_MaxN);
   sprintf(define3, "#define K_MIN %" PRIu64"\n", il_MinK);
   sprintf(define4, "#define K_MAX %" PRIu64"\n", il_MaxK);
   sprintf(define5, "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   if (ii_MaxIterations > 1)
   {
      sprintf(define6, "#define D_MAX_STEPS %d\n", ii_MaxGpuSteps);
      sprintf(define7, "#define D_MULTI_PASS\n");
   }
   
   source[0] = define1;
   source[1] = define2;
   source[2] = define3;
   source[3] = define4;
   source[4] = define5;
   
   if (ii_MaxIterations > 1)
   {
      source[5] = define6;
      source[6] = define7;
      source[7] = gfn_kernel;
      source[8] = 0;
   }
   else
   {
      source[5] = gfn_kernel;
      source[6] = 0;
   }


   ip_GFNDivisorKernel = new Kernel(ip_GFNDivisorApp->GetDevice(), "gfn_kernel", source);

   AllocatePrimeList(ip_GFNDivisorKernel->GetWorkGroupSize());

   il_RemainderList = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_FactorList  =  (uint64_t *) xmalloc(4*ii_MaxGpuFactors*sizeof(uint64_t));

   ip_KAPrime        = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   
   if (ii_MaxIterations > 1)
   {
      ip_KAParams       = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "n_range", KA_HOST_TO_GPU, (uint32_t *) &ii_Params, 5);
      ip_KARemainder    = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "remainder", KA_BIDIRECTIONAL, il_RemainderList, ii_WorkSize);
   }
   
   ip_KAFactorCount  = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList   = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);

   ip_GFNDivisorKernel->AddArgument(ip_KAPrime);
   
   if (ii_MaxIterations > 1)
   {
      ip_GFNDivisorKernel->AddArgument(ip_KAParams);
      ip_GFNDivisorKernel->AddArgument(ip_KARemainder);
   }
   
   ip_GFNDivisorKernel->AddArgument(ip_KAFactorCount);
   ip_GFNDivisorKernel->AddArgument(ip_KAFactorList);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  GFNDivisorGpuWorker::CleanUp(void)
{  
   delete ip_KAPrime;

   if (ii_MaxIterations > 1)
   {
      delete ip_KARemainder;
      delete ip_KAParams;
   }

   delete ip_KAFactorCount;
   delete ip_KAFactorList;
   delete ip_GFNDivisorKernel;

   xfree(il_RemainderList);
   xfree(il_FactorList);
}

void  GFNDivisorGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, n, iteration;
   uint64_t k, prime;
   time_t   reportTime;

   reportTime = time(NULL) + 60;

   for (iteration=0; iteration<ii_MaxIterations; iteration++)
   {
      ii_Params[0] = ii_MinN + (iteration * ii_MaxGpuSteps);
      ii_FactorCount = 0;
      
      ip_GFNDivisorKernel->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount; ii++)
      {  
         idx = ii*4;
         
         k = (uint64_t) il_FactorList[idx+0];
         n = (uint32_t) il_FactorList[idx+1];
         prime = il_FactorList[idx+2];
      
         ip_GFNDivisorApp->ReportFactor(prime, k, n, true);
         
         if (ii >= ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount, ii_MaxGpuFactors);

      if (iteration < ii_MaxIterations && ip_GFNDivisorApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_GFNDivisorApp->WriteToConsole(COT_SIEVE, "Thread %u has completed %u of %u iterations", ii_MyId, iteration, ii_MaxIterations);
         reportTime = time(NULL) + 60;
      }
   }
      
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  GFNDivisorGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GFNDivisorGpuWorker::TestMiniPrimeChunk not implemented");
}
