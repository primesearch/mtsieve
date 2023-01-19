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

   terms_t *terms = ip_MultiFactorialApp->GetTerms();
   ii_BaseCount = 0;
   
   for (idx=0; idx<ii_MultiFactorial; idx++)
      if (terms[idx].count > ii_BaseCount)
         ii_BaseCount = terms[idx].count;
   
   ii_BaseCount += 2;
   
   sprintf(defines[defineCount++], "#define D_MIN_N %u", ii_MinN);
   sprintf(defines[defineCount++], "#define D_MAX_N %u", ii_MaxN);
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define D_MULTIFACTORIAL %u", ii_MultiFactorial);
   sprintf(defines[defineCount++], "#define D_BASES %u", ii_BaseCount);
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("mf_kernel", mf_kernel, preKernelSources);

   ip_MultiFactorialApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_MultiFactorialApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   il_RemainderList = (uint64_t *) ip_Kernel->AddSharedArgument("remainders", sizeof(uint64_t), ii_PrimesInList);
   ii_Bases = (uint64_t *) ip_Kernel->AddCpuArgument("bases", sizeof(uint64_t), ii_BaseCount * ii_MultiFactorial);
   ii_Powers = (uint32_t *) ip_Kernel->AddCpuArgument("powers", sizeof(uint32_t), ii_BaseCount * ii_MultiFactorial);
   ii_Parameters = (uint32_t *) ip_Kernel->AddCpuArgument("parameters", sizeof(uint32_t), 5);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (int64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);
   
   ip_Kernel->PrintStatistics(0);

   for (idx=0; idx<ii_MultiFactorial; idx++)
   {
      // Copy this to the GPU buffer, with a fixed number of entries per idx
      memcpy(&ii_Bases[idx * ii_BaseCount], terms[idx].base, ii_BaseCount * sizeof(uint64_t));
      memcpy(&ii_Powers[idx * ii_BaseCount], terms[idx].power, ii_BaseCount * sizeof(uint32_t));
   }

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  MultiFactorialGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  MultiFactorialGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii, idx, n, firstN, lastN;
   int32_t  iteration = 0, maxIterations, c;
   uint64_t prime;
   time_t   reportTime;

   // This is an approximation.  We do not include computation of (n-1)!.
   maxIterations = ii_MultiFactorial * (ii_MaxN - ii_MinN) / ii_MaxGpuSteps;

   // For normal factorials, this loop will be iterated once.
   // For multi-factorials, it will be iterated for the value of the multi-factorial.  So
   // if this the multi-factorial is 3, i.e. n!3, then this is how it iterates:
   //    n_seed=1 will evaluate 1!3, 4!3, 7!3, etc.  (7!3 = 7*4*1)
   //    n_seed=2 will evaluate 2!3, 5!3, 8!3, etc.  (8!3 = 8*5*2)
   //    n_seed=3 will evaluate 3!3, 6!3, 9!3, etc.  (9!3 = 9*6*3)

   for (uint32_t mf=0; mf<ii_MultiFactorial; mf++)
   {
      // If ii_Multifactorial is even and mf is odd then 
      // n!ii_Multifactorial+1 and n!ii_Multifactorial-1 are always even
      // when mf is odd, so we do not need to go any further.
      if (!(ii_MultiFactorial & 1) && (mf & 1))
         continue;

      reportTime = time(NULL) + 60;

      ii_Parameters[0] = mf;
      ii_Parameters[1] = 1;
      
      firstN = ii_MinN - 1;
      
      while (firstN % ii_MultiFactorial != mf)
         firstN--;

      do
      {
         iteration++;
         lastN = firstN + (ii_MaxGpuSteps * ii_MultiFactorial);

         ii_Parameters[2] = firstN;
         ii_Parameters[3] = lastN;
         
         if (ii_Parameters[3] > ii_MaxN)
            ii_Parameters[3] = ii_MaxN;

         ii_FactorCount[0] = 0;
        
         ip_Kernel->Execute(ii_PrimesInList);

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
      
         firstN = lastN;
         ii_Parameters[1] = 0;
      } while (firstN < ii_MaxN);
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  MultiFactorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("MultiFactorialGpuWorker::TestMiniPrimeChunk not implemented");
}
