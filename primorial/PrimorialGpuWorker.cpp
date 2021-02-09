/* PrimorialGpuWorker.cpp -- (C) Mark Rodenkirch, September 2016

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "PrimorialGpuWorker.h"
#include "primorial_kernel.h"
#include <time.h>
#include "../sieve/primesieve.hpp"

PrimorialGpuWorker::PrimorialGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        define1[50];
   char        define2[50];
   char        define3[50];
   char        define4[50];
   char        define5[50];
   char        define6[50];
   char        define7[50];
   const char *source[10];
   uint32_t    idx, minIdx = 0;
   uint16_t    biggestGap = 0;

   vector<uint64_t>  primes;
   vector<uint64_t>::iterator it;
   
   ib_GpuWorker = true;
   
   ip_PrimorialApp = (PrimorialApp *) theApp;

   ii_MinPrimorial = ip_PrimorialApp->GetMinPrimorial();
   ii_MaxPrimorial = ip_PrimorialApp->GetMaxPrimorial();
   
   ii_MaxGpuSteps = ip_PrimorialApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_PrimorialApp->GetMaxGpuFactors();
   
   ip_PrimorialPrimes = ip_PrimorialApp->GetPrimorialPrimes(ii_NumberOfPrimorialPrimes);
   ip_PrimorialPrimeGaps = ip_PrimorialApp->GetPrimorialPrimeGaps(biggestGap);
   
   for (idx=0; idx<ii_NumberOfPrimorialPrimes; idx++)
      if (ii_MinPrimorial == ip_PrimorialPrimes[idx])
         minIdx = idx;
   
   ii_MaxIterations = 1 + ((1 + ii_NumberOfPrimorialPrimes) / ii_MaxGpuSteps);

   sprintf(define1, "#define D_MIN_IDX %d\n", minIdx);
   sprintf(define2, "#define D_MAX_STEPS %d\n", ii_MaxGpuSteps);
   sprintf(define3, "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   sprintf(define4, "#define D_BIGGEST_GAP %u\n", biggestGap);
   sprintf(define5, "#define D_MIN_PRIMORIAL %u\n", ii_MinPrimorial);
   sprintf(define6, "#define FIRST_PRIMORIAL %u\n", FIRST_PRIMORIAL);
   sprintf(define7, "#define FIRST_PRIMORIAL_PRIME %u\n", FIRST_PRIMORIAL_PRIME);
   
   source[0] = define1;
   source[1] = define2;
   source[2] = define3;
   source[3] = define4;
   source[4] = define5;
   source[5] = define6;
   source[6] = define7;
   source[7] = primorial_kernel;
   source[8] = 0;

   ip_PrimorialKernel = new Kernel(ip_PrimorialApp->GetDevice(), "primorial_kernel", source);

   AllocatePrimeList(ip_PrimorialKernel->GetWorkGroupSize());
   
   ip_RemainderList = (uint64_t *) xmalloc(2*ii_WorkSize*sizeof(uint64_t));
   ip_FactorList    = (int64_t *)  xmalloc(4*ii_MaxGpuFactors*sizeof(int64_t));

   ip_KAPrime         = new KernelArgument(ip_PrimorialApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KAPrimorialGaps = new KernelArgument(ip_PrimorialApp->GetDevice(), "primeGaps", KA_HOST_TO_GPU, ip_PrimorialPrimeGaps, ii_NumberOfPrimorialPrimes);
   ip_KAParams        = new KernelArgument(ip_PrimorialApp->GetDevice(), "n_range", KA_HOST_TO_GPU, (uint64_t *) &ii_Params, 5);
   ip_KARemainder     = new KernelArgument(ip_PrimorialApp->GetDevice(), "remainder", KA_BIDIRECTIONAL, ip_RemainderList, 2*ii_WorkSize);
   ip_KAFactorCount   = new KernelArgument(ip_PrimorialApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList    = new KernelArgument(ip_PrimorialApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, ip_FactorList, 4*ii_MaxGpuFactors);

   ip_PrimorialKernel->AddArgument(ip_KAPrime);
   ip_PrimorialKernel->AddArgument(ip_KAPrimorialGaps);
   ip_PrimorialKernel->AddArgument(ip_KAParams);
   ip_PrimorialKernel->AddArgument(ip_KARemainder);
   ip_PrimorialKernel->AddArgument(ip_KAFactorCount);
   ip_PrimorialKernel->AddArgument(ip_KAFactorList);
   
   ip_PrimorialKernel->PrintStatistics(0);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  PrimorialGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KAPrimorialGaps;
   delete ip_KARemainder;
   delete ip_KAParams;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;

   xfree(il_PrimeList);
   xfree(ip_RemainderList);
   xfree(ip_FactorList);
}

void  PrimorialGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii, idx, primorial;
   uint32_t iteration;
   int32_t  c;
   uint64_t primeFactor;
   time_t   reportTime;

   reportTime = time(NULL) + 60;

   for (iteration=0; iteration<ii_MaxIterations; iteration++)
   {
      // Set the starting index and primorial for this iteration
      ii_Params[0] = iteration * ii_MaxGpuSteps;

      ii_FactorCount = 0;
  
      ip_PrimorialKernel->Execute(ii_WorkSize);
      
      for (ii=0; ii<ii_FactorCount; ii++)
      {  
         idx = ii*4;
         
         primorial = (uint32_t) ip_FactorList[idx+0];
         c = (int32_t) ip_FactorList[idx+1];
         primeFactor = ip_FactorList[idx+2];
      
         ip_PrimorialApp->ReportFactor(primeFactor, primorial, c, true);
            
         if (ii >= ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount, ii_MaxGpuFactors);
            
      if (iteration < ii_MaxIterations && ip_PrimorialApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_PrimorialApp->WriteToConsole(COT_SIEVE, "Thread %u has completed %u of %u iterations", ii_MyId, iteration, ii_MaxIterations);
         reportTime = time(NULL) + 60;
      }
   }

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  PrimorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("PrimorialGpuWorker::TestMiniPrimeChunk not implemented");
}
