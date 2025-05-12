/* HyperCullenWoodallGpuWorker.cpp -- (C) Mark Rodenkirch

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "HyperCullenWoodallGpuWorker.h"
#include "hcw_kernel.gpu.h"

HyperCullenWoodallGpuWorker::HyperCullenWoodallGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ib_GpuWorker = true;
   
   ip_HyperCullenWoodallApp = (HyperCullenWoodallApp *) theApp;

   // Allocate enough memory to hold all of the terms.
   ii_GroupSize = ip_HyperCullenWoodallApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_HyperCullenWoodallApp->GetMaxGpuFactors();

   if (ip_HyperCullenWoodallApp->IsPlus())
      snprintf(defines[defineCount++], 50, "#define IS_PLUS");
   
   if (ip_HyperCullenWoodallApp->IsMinus())
      snprintf(defines[defineCount++], 50, "#define IS_MINUS");

   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);
      
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("hcw_kernel", hcw_kernel, preKernelSources);

   ip_HyperCullenWoodallApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_HyperCullenWoodallApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   ii_KernelTerms = (uint32_t *) ip_Kernel->AddCpuArgument("terms", sizeof(uint32_t), ii_GroupSize);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   // The thread can't start until initialization is done
   ib_Initialized = true;

   il_NextTermsBuild = 0;
   ip_GpuGroup = NULL;
}

void  HyperCullenWoodallGpuWorker::CleanUp(void)
{
   ip_HyperCullenWoodallApp->FreeGpuGroupsOfTerms(ip_GpuGroup);

   delete ip_Kernel;
}

void  HyperCullenWoodallGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t    group, groups;
   uint32_t    idx;
   uint32_t    x, y;
   int32_t     sign;
   uint64_t    thePrime;
   time_t      reportTime = time(NULL) + 60;
   gpugroup_t *gPtr;

   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      ip_HyperCullenWoodallApp->FreeGpuGroupsOfTerms(ip_GpuGroup);
      
      ip_GpuGroup = ip_HyperCullenWoodallApp->GetGpuGroupsOfTerms();
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
   }

   groups = 0;
   gPtr = ip_GpuGroup;
   while (gPtr != NULL)
   {
      groups++;
      gPtr = (gpugroup_t *) gPtr->nextGroup;
   }

   gPtr = ip_GpuGroup;
   group = 0;
   while (gPtr != NULL)
   {
      if (ip_HyperCullenWoodallApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_HyperCullenWoodallApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d pf %d iterations", ii_MyId, group, groups);
         reportTime = time(NULL) + 60;
      }
            
      memcpy(ii_KernelTerms, gPtr->terms, ii_GroupSize * sizeof(uint32_t));
      
      ii_FactorCount[0] = 0;
      
      ip_Kernel->Execute(ii_PrimesInList);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {  
         idx = ii*4;
         
         x = (uint32_t) il_FactorList[idx+0];
         y = (uint32_t) il_FactorList[idx+1];
         sign = (il_FactorList[idx+2] == +1 ? +1 : -1);
         thePrime = il_FactorList[idx+3];

         ip_HyperCullenWoodallApp->ReportFactor(thePrime, x, y, sign);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);

      group++;
      gPtr = (gpugroup_t *) gPtr->nextGroup;      
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  HyperCullenWoodallGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("HyperCullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}
