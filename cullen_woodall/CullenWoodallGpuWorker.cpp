/* CullenWoodallGpuWorker.cpp -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CullenWoodallGpuWorker.h"
#include "cw_kernel.gpu.h"

CullenWoodallGpuWorker::CullenWoodallGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx, groupCount;
   
   ib_GpuWorker = true;
   
   ip_CullenWoodallApp = (CullenWoodallApp *) theApp;

   ii_Base = ip_CullenWoodallApp->GetBase();
   ii_MaxGpuSteps = ip_CullenWoodallApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_CullenWoodallApp->GetMaxGpuFactors();

   groupCount = 1 + (ip_CullenWoodallApp->GetTermCount() / ii_MaxGpuSteps);

   // Allocate enough memory to hold all of the terms.
   ii_Terms = (uint32_t *) xmalloc(groupCount*ii_MaxGpuSteps*sizeof(int32_t));
   
   sprintf(defines[defineCount++], "#define BASE %d", ii_Base);
   
   if (ip_CullenWoodallApp->IsCullenSearch())
      sprintf(defines[defineCount++], "#define CHECK_CULLEN");
   
   if (ip_CullenWoodallApp->IsWoodallSearch())
      sprintf(defines[defineCount++], "#define CHECK_WOODALL");
   
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %d", ii_MaxGpuFactors);
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("cw_kernel", cw_kernel, preKernelSources);

   ip_CullenWoodallApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_CullenWoodallApp->GetGpuPrimesPerWorker();

   // Add space for a few extra terms to guarantee that the last term in the group is 0.
   ii_GroupSize = ii_MaxGpuSteps + 10;
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   ii_KernelTerms = (uint32_t *) ip_Kernel->AddCpuArgument("terms", sizeof(uint32_t), ii_GroupSize);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (int64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(0);
   
   il_NextTermsBuild = 0;
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CullenWoodallGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  CullenWoodallGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t group, idx;
   int32_t  n, c;
   uint64_t prime;
   
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      ii_Groups = ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MaxGpuSteps, ii_GroupSize);
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
   }

   for (group=0; group<ii_Groups; group++)
   {
      ii_FactorCount[0] = 0;
      idx = group * ii_GroupSize;
      
      memcpy(ii_KernelTerms, &ii_Terms[idx], ii_GroupSize * sizeof(uint32_t));

      ip_Kernel->Execute(ii_PrimesInList);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {  
         idx = ii*4;
         
         n = (uint32_t) il_FactorList[idx+0];
         c = (int32_t) il_FactorList[idx+1];
         prime = il_FactorList[idx+2];

         ip_CullenWoodallApp->ReportFactor(prime, n, c);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
   }

   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  CullenWoodallGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}
