/* AlternatingFactorial.cpp -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "AlternatingFactorialGpuWorker.h"
#include "af_kernel.gpu.h"

AlternatingFactorialGpuWorker::AlternatingFactorialGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    idx, defineCount;
   
   ib_GpuWorker = true;
   
   ip_AlternatingFactorialApp = (AlternatingFactorialApp *) theApp;
   
   ii_MaxGpuFactors = ip_AlternatingFactorialApp->GetMaxGpuFactors() + 10;
   
   defineCount = 0;
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define D_MAX_N %u", ip_AlternatingFactorialApp->GetMaxN());
   sprintf(defines[defineCount++], "#define D_MAX_STEPS %u", ip_AlternatingFactorialApp->GetMaxGpuSteps());
      
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("af_kernel", af_kernel, preKernelSources);

   ip_AlternatingFactorialApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_AlternatingFactorialApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   ii_Parameters = (uint32_t *) ip_Kernel->AddCpuArgument("parameters", sizeof(uint32_t), 4);
   il_FactorialResiduals = (uint64_t *) ip_Kernel->AddSharedArgument("factorialResiduals", sizeof(uint64_t), ii_WorkSize);
   il_AltFactorialResiduals = (uint64_t *) ip_Kernel->AddSharedArgument("altFactorialResiduals", sizeof(uint64_t), ii_WorkSize);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 2*ii_MaxGpuFactors);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  AlternatingFactorialGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  AlternatingFactorialGpuWorker::TestMegaPrimeChunk(void)
{
   uint64_t thePrime;
   uint32_t term, idx;

   ii_Parameters[0] = 1;

   do
   {      
      ii_FactorCount[0] = 0;
      
      ip_Kernel->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {
         idx = ii*2;
         
         term = (uint32_t) il_FactorList[idx+0];
         thePrime = il_FactorList[idx+1];
         
         ip_AlternatingFactorialApp->ReportFactor(thePrime, term);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }
   
      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
      
      ii_Parameters[0] += ip_AlternatingFactorialApp->GetMaxGpuSteps();
   } while (ii_Parameters[0] < ip_AlternatingFactorialApp->GetMaxN());

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  AlternatingFactorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("AlternatingFactorialGpuWorker::TestMiniPrimeChunk not implemented");
}
