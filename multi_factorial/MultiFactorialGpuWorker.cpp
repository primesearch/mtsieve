/* MultiFactorialGpuWorker.cpp -- (C) Mark Rodenkirch, September 2016

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include <cinttypes>
#include "MultiFactorialGpuWorker.h"
#include "mf_kernel.gpu.h"

MultiFactorialGpuWorker::MultiFactorialGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;

   ib_GpuWorker = true;
   
   ip_MultiFactorialApp = (MultiFactorialApp *) theApp;

   ii_MinN = ip_MultiFactorialApp->GetMinN();
   ii_MaxN = ip_MultiFactorialApp->GetMaxN();
   ii_MaxGpuSteps = ip_MultiFactorialApp->GetMaxGpuSteps();
   ii_MultiFactorial = ip_MultiFactorialApp->GetMultiFactorial();
   
   ii_MaxGpuFactors = ip_MultiFactorialApp->GetMaxGpuFactors();

   sprintf(defines[defineCount++], "#define D_MIN_N %d", ii_MinN);
   sprintf(defines[defineCount++], "#define D_MAX_N %d", ii_MaxN);
   sprintf(defines[defineCount++], "#define D_MAX_STEPS %d", ii_MaxGpuSteps);
   sprintf(defines[defineCount++], "#define D_MULTIFACTORIAL %d", ii_MultiFactorial);
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %d", ii_MaxGpuFactors);
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("mf_kernel", mf_kernel, preKernelSources);

   ip_MultiFactorialApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_MultiFactorialApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   il_RemainderList = (uint64_t *) ip_Kernel->AddGpuArgument("remainders", sizeof(uint64_t), 2*ii_WorkSize);
   ii_Parameters = (uint32_t *) ip_Kernel->AddCpuArgument("parameters", sizeof(uint32_t), 5);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (int64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);
   
   ip_Kernel->PrintStatistics(0);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  MultiFactorialGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  MultiFactorialGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii, idx;
   uint32_t startN, n;
   int32_t  iteration = 0, maxIterations, c;
   uint64_t prime;
   time_t   reportTime;

   // This is an approximation
   maxIterations = ii_MultiFactorial * (1 + ii_MaxN) / ii_MaxGpuSteps;

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
   for (startN=1; startN<=ii_MultiFactorial; startN++)
   {
      // If startN is odd and mf is even, then i!mf is always odd, thus
      // startN!mf+1 and startN!mf-1 are always even.
      if (!(ii_MultiFactorial & 1) && (startN & 1))
         continue;

      // The first parameter is the starting n for the calculation, i.e. n!
      ii_Parameters[0] = startN;
      ii_Parameters[1] = startN;

      reportTime = time(NULL) + 60;
      
      do
      {
         iteration++;
         ii_FactorCount[0] = 0;
     
         ip_Kernel->Execute(ii_WorkSize);

         for (ii=0; ii<ii_FactorCount[0]; ii++)
         {  
            idx = ii*4;
            
            n = (uint32_t) il_FactorList[idx+0];
            c = (int32_t) il_FactorList[idx+1];
            prime = il_FactorList[idx+2];
         
            ip_MultiFactorialApp->ReportFactor(prime, n, c);
            
            if ((ii+1) == ii_MaxGpuFactors)
               break;
         }

         if (ii_FactorCount[0] >= ii_MaxGpuFactors)
            FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);

         if (iteration < maxIterations && ip_MultiFactorialApp->IsInterrupted() && time(NULL) > reportTime)
         {
            ip_MultiFactorialApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, iteration, maxIterations);
            reportTime = time(NULL) + 60;
         }
         
         // Set where the next range is starting.
         ii_Parameters[1] += (ii_MaxGpuSteps * ii_MultiFactorial);
      } while (ii_Parameters[1] <= ii_MaxN);
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  MultiFactorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("MultiFactorialGpuWorker::TestMiniPrimeChunk not implemented");
}
