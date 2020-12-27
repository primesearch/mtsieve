/* CullenWoodallGpuWorker.cpp -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CullenWoodallGpuWorker.h"
#include "cw_kernel.h"
#include "../x86_asm/fpu-asm-x86.h"

CullenWoodallGpuWorker::CullenWoodallGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   const char *gcwSource[9];
   char  maxFactors[50];
   char  isCullen[50];
   char  isWoodall[50];
   char  theBase[50];
   
   ib_GpuWorker = true;
   
   ip_CullenWoodallApp = (CullenWoodallApp *) theApp;

   ii_Base = ip_CullenWoodallApp->GetBase();
   ii_MaxGpuSteps = ip_CullenWoodallApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_CullenWoodallApp->GetMaxGpuFactors();

   sprintf(theBase, "#define BASE %d\n", ii_Base);
   
   if (ip_CullenWoodallApp->IsCullenSearch())
      sprintf(isCullen, "#define CHECK_CULLEN\n"); 
   else
      sprintf(isCullen, "\n");
   
   if (ip_CullenWoodallApp->IsWoodallSearch())
      sprintf(isWoodall, "#define CHECK_WOODALL\n");
   else
      sprintf(isWoodall, "\n");
   
   sprintf(maxFactors, "#define MAX_FACTORS %d\n", ii_MaxGpuFactors);

   gcwSource[0] = theBase;
   gcwSource[1] = isCullen;
   gcwSource[2] = isWoodall;
   gcwSource[3] = maxFactors;
   gcwSource[4] = cw_kernel;
   gcwSource[5] = 0;

   ip_GCWKernel = new Kernel(ip_CullenWoodallApp->GetDevice(), "cw_kernel", gcwSource);
   
   AllocatePrimeList(ip_GCWKernel->GetWorkGroupSize());
   
   // If there are more y for any x than the number of steps, then we'll give an error
   // to the user at runtime when building the term groups.
   ii_GroupSize = ii_MaxGpuSteps + 10;

   uint32_t groupCount = 1 + (ip_CullenWoodallApp->GetTermCount() / ii_MaxGpuSteps);
   
   // Allocate enough memory to hold all of the terms.
   ii_Terms = (uint32_t *) xmalloc(groupCount*ii_MaxGpuSteps*sizeof(int32_t));

   il_FactorList  =  (int64_t *) xmalloc(4*ii_MaxGpuFactors*sizeof(int64_t));

   ip_KAPrime        = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);

   // We will change the pointer before executing the kernel
   ip_KATerms          = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "terms", KA_HOST_TO_GPU, ii_Terms, ii_GroupSize);
   ip_KAFactorCount    = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList     = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);

   ip_GCWKernel->AddArgument(ip_KAPrime);
   ip_GCWKernel->AddArgument(ip_KATerms);
   ip_GCWKernel->AddArgument(ip_KAFactorCount);
   ip_GCWKernel->AddArgument(ip_KAFactorList);

   // The thread can't start until initialization is done
   ib_Initialized = true;

   il_NextTermsBuild = 0;
}

void  CullenWoodallGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;

   delete ip_KATerms;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;

   delete ip_GCWKernel;
   
   xfree(il_FactorList);
   xfree(ii_Terms);
}

void  CullenWoodallGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t group, idx;
   int32_t  n, c;
   uint64_t prime;
   time_t   reportTime = time(NULL) + 60;
   
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      ii_Groups = ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MaxGpuSteps, ii_GroupSize);
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
   }
      
   for (group=0; group<ii_Groups; group++)
   {
      if (ip_CullenWoodallApp->IsInterrupted() || time(NULL) > reportTime)
      {
         ip_CullenWoodallApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, group, ii_Groups);
         reportTime = time(NULL) + 60;
      }

      ii_FactorCount = 0;
      idx = group * ii_GroupSize;
      ip_KATerms->SetHostMemory(&ii_Terms[idx]);

      ip_GCWKernel->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount; ii++)
      {  
         idx = ii*4;
         
         n = (uint32_t) il_FactorList[idx+0];
         c = (int32_t) il_FactorList[idx+1];
         prime = il_FactorList[idx+2];
      
         if (ip_CullenWoodallApp->ReportFactor(prime, n, c))
            VerifyFactor(prime, n, c);
      }

      if (ii_FactorCount >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors.  Use -M to increase max factors", ii_FactorCount);
   }

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  CullenWoodallGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}

void  CullenWoodallGpuWorker::VerifyFactor(uint64_t p, uint32_t n, int32_t c)
{
   uint64_t rem;
   
   fpu_push_1divp(p);
      
   rem = fpu_powmod(ip_CullenWoodallApp->GetBase(), n, p);
   rem = fpu_mulmod(rem, n, p);
   
   fpu_pop();
   
   if (c == -1 && rem != +1)
      FatalError("%" PRIu64" does not divide %u*%u^%u-1", p, n, ii_Base, n);
   
   if (c == +1 && rem != p-1)
      FatalError("%" PRIu64" does not divide %u*%u^%u+1", p, n, ii_Base, n);
}
