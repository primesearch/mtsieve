/* SmarandacheGpuWorker.cpp -- (C) Mark Rodenkirch, January 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "SmarandacheGpuWorker.h"
#include "mf_kernel.h"
#include "../x86_asm/fpu-asm-x86.h"
#include <time.h>

SmarandacheGpuWorker::SmarandacheGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        define1[40];
   char        define2[40];
   char        define3[40];
   char        define4[40];
   char        define5[40];
   const char *source[10];

   ib_GpuWorker = true;
   
   ip_SmarandacheApp = (SmarandacheApp *) theApp;

   ii_MinN = ip_SmarandacheApp->GetMinN();
   ii_MaxN = ip_SmarandacheApp->GetMaxN();
   ii_MaxGpuSteps = ip_SmarandacheApp->GetMaxGpuSteps();
   ii_Smarandache = ip_SmarandacheApp->GetSmarandache();
   
   ii_MaxGpuFactors = ip_SmarandacheApp->GetMaxGpuFactors();
   
   sprintf(define1, "#define D_MIN_N %d\n", ii_MinN);
   sprintf(define2, "#define D_MAX_N %d\n", ii_MaxN);
   sprintf(define3, "#define D_MAX_STEPS %d\n", ii_MaxGpuSteps);
   sprintf(define4, "#define D_Smarandache %d\n", ii_Smarandache);
   sprintf(define5, "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   
   source[0] = define1;
   source[1] = define2;
   source[2] = define3;
   source[3] = define4;
   source[4] = define5;
   source[5] = mf_kernel;
   source[6] = 0;

   ip_FactorialKernel = new Kernel(ip_SmarandacheApp->GetDevice(), "mf_kernel", source);

   AllocatePrimeList(ip_FactorialKernel->GetWorkGroupSize());
   
   il_RemainderList = (uint64_t *) xmalloc(2*ii_WorkSize*sizeof(uint64_t));
   il_FactorList    = (int64_t *)  xmalloc(4*ii_MaxGpuFactors*sizeof(int64_t));

   ip_KAPrime        = new KernelArgument(ip_SmarandacheApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KARemainder    = new KernelArgument(ip_SmarandacheApp->GetDevice(), "remainder", KA_BIDIRECTIONAL, il_RemainderList, 2*ii_WorkSize);
   ip_KAParams       = new KernelArgument(ip_SmarandacheApp->GetDevice(), "n_range", KA_HOST_TO_GPU, (int32_t *) &ii_Params, 5);
   ip_KAFactorCount  = new KernelArgument(ip_SmarandacheApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList   = new KernelArgument(ip_SmarandacheApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);

   ip_FactorialKernel->AddArgument(ip_KAPrime);
   ip_FactorialKernel->AddArgument(ip_KARemainder);
   ip_FactorialKernel->AddArgument(ip_KAParams);
   ip_FactorialKernel->AddArgument(ip_KAFactorCount);
   ip_FactorialKernel->AddArgument(ip_KAFactorList);
   
   ip_FactorialKernel->PrintStatistics(0);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  SmarandacheGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KARemainder;
   delete ip_KAParams;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;

   xfree(il_PrimeList);
   xfree(il_RemainderList);
   xfree(il_FactorList);
}

void  SmarandacheGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii, idx;
   uint32_t startN, n;
   int32_t  iteration = 0, maxIterations, c;
   uint64_t prime;
   time_t   reportTime;

   // This is an approximation
   maxIterations = ii_Smarandache * (1 + ii_MaxN) / ii_MaxGpuSteps;

   // For normal factorials, this loop will be iterated one.
   // For multi-factorials, it will be iterated for the value of the multi-factorial.  So
   // if this the multi-factorial is 3, i.e. n!3, then this is how it iterates:
   //    n_seed=1 will evaluate 1!3, 4!3, 7!3, etc.  (7!3 = 7*4*1)
   //    n_seed=2 will evaluate 2!3, 5!3, 8!3, etc.  (8!3 = 8*5*2)
   //    n_seed=3 will evaluate 3!3, 6!3, 9!3, etc.  (9!3 = 9*6*3)
   // For inverse multi-factorials, it will be similar.  If the inverse multi-factorial
   // is 3, i.e. n!/n!3, then this is how it iterates:
   //    n_seed=1 will evaluate 1!/1!3, 4!/4!3, 7!/7!3, etc.  (7!/7!3 = 6*5*3*2, i.e not multiplying by 7, 4, and 1)
   //    n_seed=2 will evaluate 2!/2!3, 5!/5!3, 8!/8!3, etc.  (8!/8!3 = 7*6*4*3*1, i.e. not multiplying by 8, 5, and 2)
   //    n_seed=3 will evaluate 3!/3!3, 6!/6!3, 9!/9!3, etc.  (9!/9!3 = 8*7*5*4*2*1, i.e. not multiplying by 9, 6, and 3))
   for (startN=1; startN<=ii_Smarandache; startN++)
   {
      // If startN is odd and mf is even, then i!mf is always odd, thus
      // startN!mf+1 and startN!mf-1 are always even.
      if (!(ii_Smarandache & 1) && (startN & 1))
         continue;

      // The first parameter is the starting n for the calculation, i.e. n!
      ii_Params[0] = startN;

      // Initialize the list of remainders
      for (idx=0; idx<ii_WorkSize; idx++)
         il_RemainderList[idx] = 1;

      reportTime = time(NULL) + 60;

      ii_Params[1] = startN;
      
      do
      {
         iteration++;
         ii_FactorCount = 0;
     
         ip_FactorialKernel->Execute(ii_WorkSize);

         for (ii=0; ii<ii_FactorCount; ii++)
         {  
            idx = ii*4;
            
            n = (uint32_t) il_FactorList[idx+0];
            c = (int32_t) il_FactorList[idx+1];
            prime = il_FactorList[idx+2];
         
            ip_SmarandacheApp->ReportFactor(prime, n, c);
            
            if (ii >= ii_MaxGpuFactors)
               break;
         }

         if (ii_FactorCount >= ii_MaxGpuFactors)
            FatalError("Could not handle all GPU factors.  A range of p generated %u factors.  Use -M to increase max factors", ii_FactorCount);

         if (iteration < maxIterations && ip_SmarandacheApp->IsInterrupted() && time(NULL) > reportTime)
         {
            ip_SmarandacheApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, iteration, maxIterations);
            reportTime = time(NULL) + 60;
         }
         
         // Set where the next range is starting.
         ii_Params[1] += (ii_MaxGpuSteps * ii_Smarandache);
      } while (ii_Params[1] <= ii_MaxN);
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  SmarandacheGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheGpuWorker::TestMiniPrimeChunk not implemented");
}
