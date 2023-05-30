/* XYYXSparseGpuWorker.cpp -- (C) Mark Rodenkirch, September 2012

   This class sets up the call to the XYYXSparseGpuWorker GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "XYYXSparseGpuWorker.h"
#include "xyyx_sparse_kernel.gpu.h"

XYYXSparseGpuWorker::XYYXSparseGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ib_GpuWorker = true;
   
   ip_XYYXApp = (XYYXApp *) theApp;

   ii_GroupSize = ip_XYYXApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_XYYXApp->GetMaxGpuFactors();
   ip_Terms = 0;
   
   if (ip_XYYXApp->IsPlus())
      snprintf(defines[defineCount++], 50, "#define IS_PLUS");
   
   if (ip_XYYXApp->IsMinus())
      snprintf(defines[defineCount++], 50, "#define IS_MINUS");

   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("xyyx_sparse_kernel", xyyx_sparse_kernel, preKernelSources);

   ip_XYYXApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_XYYXApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   ii_KernelTerms = (uint32_t *) ip_Kernel->AddCpuArgument("terms", sizeof(uint32_t), 2*ii_GroupSize);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   // The thread can't start until initialization is done
   ib_Initialized = true;

   il_NextTermsBuild = 0;
}

void  XYYXSparseGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
   
   xfree(ip_Terms);
}

void  XYYXSparseGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t  group;
   uint32_t  idx;
   uint32_t  x, y;
   uint64_t  thePrime;
   time_t    reportTime = time(NULL) + 60;

   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      if (il_NextTermsBuild > 0)
         xfree(ip_Terms);
      
      ip_Terms = ip_XYYXApp->GetSparseGroupedTerms();
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
      
      ii_Groups = 0;
      
      for (uint32_t tIdx=0; ; tIdx+=ii_GroupSize)
      {
         if (ip_Terms[tIdx].x == 0)
            break;
         
         ii_Groups++;
      }      
   }

   for (group=0; group<ii_Groups; group++)
   {
      if (ip_XYYXApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_XYYXApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, group, ii_Groups);
         reportTime = time(NULL) + 60;
      }
      
      idx = group * ii_GroupSize;
      
      memcpy(ii_KernelTerms, &ip_Terms[idx], ii_GroupSize * sizeof(gputerm_t));
      
      ii_FactorCount[0] = 0;
      
      ip_Kernel->Execute(ii_PrimesInList);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {  
         idx = ii*4;
         
         x = (uint32_t) il_FactorList[idx+0];
         y = (uint32_t) il_FactorList[idx+1];
         thePrime = il_FactorList[idx+2];

         ip_XYYXApp->ReportFactor(thePrime, x, y);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  XYYXSparseGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("XYYXSparseGpuWorker::TestMiniPrimeChunk not implemented");
}
