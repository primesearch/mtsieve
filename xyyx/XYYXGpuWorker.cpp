/* XYYXGpuWorker.cpp -- (C) Mark Rodenkirch, September 2012

   This class sets up the call to the XYYXGpuWorker GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "XYYXGpuWorker.h"
#include "xyyx_kernel.gpu.h"

XYYXGpuWorker::XYYXGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ib_GpuWorker = true;
   
   ip_XYYXApp = (XYYXApp *) theApp;

   // Allocate enough memory to hold all of the terms.
   ii_Groups = ip_XYYXApp->GetNumberOfGroups();
   ii_GroupSize = ip_XYYXApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_XYYXApp->GetMaxGpuFactors();

   ii_Terms = (uint32_t *) xmalloc(ii_Groups * ii_GroupSize * sizeof(uint32_t));
   
   if (ip_XYYXApp->IsPlus())
      sprintf(defines[defineCount++], "#define IS_PLUS");
   
   if (ip_XYYXApp->IsMinus())
      sprintf(defines[defineCount++], "#define IS_MINUS");

   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);
      
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = new Kernel(ip_XYYXApp->GetDevice(), "xyyx_kernel", xyyx_kernel, preKernelSources);

   ip_XYYXApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_XYYXApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   ii_KernelTerms = (uint32_t *) ip_Kernel->AddCpuArgument("terms", sizeof(uint32_t), ii_GroupSize);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   // The thread can't start until initialization is done
   ib_Initialized = true;

   il_NextTermsBuild = 0;
}

void  XYYXGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  XYYXGpuWorker::TestMegaPrimeChunk(void)
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
      ii_Groups = ip_XYYXApp->GetGroupedTerms(ii_Terms);
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
   }

   for (group=0; group<ii_Groups; group++)
   {
      if (ip_XYYXApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_XYYXApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, group, ii_Groups);
         reportTime = time(NULL) + 60;
      }
      
      idx = group * ii_GroupSize;
      
      memcpy(ii_KernelTerms, &ii_Terms[idx], ii_GroupSize * sizeof(uint32_t));
      
      ii_FactorCount[0] = 0;
      
      ip_Kernel->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {  
         idx = ii*4;
         
         x = (uint32_t) il_FactorList[idx+0];
         y = (int32_t) il_FactorList[idx+1];
         thePrime = il_FactorList[idx+2];

         ip_XYYXApp->ReportFactor(thePrime, x, y);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  XYYXGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("XYYXGpuWorker::TestMiniPrimeChunk not implemented");
}
