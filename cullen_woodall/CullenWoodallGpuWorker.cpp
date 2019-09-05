/* CullenWoodallGpuWorker.cpp -- (C) Mark Rodenkirch, May 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CullenWoodallGpuWorker.h"
#include "magiccl.h"
#include "cw_kernel.h"
#include "../x86_asm/fpu-asm-x86.h"

CullenWoodallGpuWorker::CullenWoodallGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ib_GpuWorker = true;
   
   ip_CullenWoodallApp = (CullenWoodallApp *) theApp;

   ii_Base = ip_CullenWoodallApp->GetBase();

   CreateKernels(false);

   AllocatePrimeList(ip_GCWKernel->GetWorkGroupSize());
      
   // Assume that fewer than 1 in 100 primes will factor one of the remaining terms.
   // Since the first chunk of primes is sieved only with the CPU, this means that
   // the remaining terms for the first GPU prime chunk should reveal a very small
   // percentage of factors.
   ii_MaxFactorCount = ii_WorkSize * 2;
   
   // Delete then re-create the kernels with the worksize specified in #defines
   delete(ip_MagicKernel);
   delete(ip_GCWKernel);
   
   CreateKernels(true);
   
   // If there are more y for any x than the number of steps, then we'll give an error
   // to the user at runtime when building the term groups.
   ii_MinGroupSize = ip_CullenWoodallApp->GetGpuSteps();
   ii_MaxGroupSize = 2 * ii_MinGroupSize;

   ii_MaxTermCount = ip_CullenWoodallApp->GetMaxN() - ip_CullenWoodallApp->GetMinN() + 10;
   
   // Ensure the buffer big enough to avoid issues when executing the kernel
   if (ii_MaxTermCount < ii_MaxGroupSize)
      ii_MaxTermCount = ii_MaxGroupSize;

   // Allocate enough memory to hold all of the terms.
   ii_Terms = (uint32_t *) xmalloc(ii_MaxTermCount*sizeof(int32_t));

   il_MagicNumber = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MagicShift  = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_FactorList  =  (int64_t *) xmalloc(4*ii_MaxFactorCount*sizeof(int64_t));

   ip_KAPrime        = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   
   ip_MKAMagicNumber = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "magic_number", KA_GPU_TO_HOST, il_MagicNumber, ii_WorkSize);
   ip_MKAMagicShift  = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "magic_shift", KA_GPU_TO_HOST, il_MagicShift, ii_WorkSize);

   ip_GCWKAMagicNumber = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "magic_number", KA_HOST_TO_GPU, il_MagicNumber, ii_WorkSize);
   ip_GCWKAMagicShift  = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "magic_shit", KA_HOST_TO_GPU, il_MagicShift, ii_WorkSize);
   
   // We will change the pointer before executing the kernel
   ip_KATerms          = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "terms", KA_HOST_TO_GPU, ii_Terms, ii_MaxGroupSize);
   ip_KAFactorCount    = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList     = new KernelArgument(ip_CullenWoodallApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxFactorCount);

   ip_MagicKernel->AddArgument(ip_KAPrime);
   ip_MagicKernel->AddArgument(ip_MKAMagicNumber);
   ip_MagicKernel->AddArgument(ip_MKAMagicShift);

   ip_GCWKernel->AddArgument(ip_KAPrime);
   ip_GCWKernel->AddArgument(ip_GCWKAMagicNumber);
   ip_GCWKernel->AddArgument(ip_GCWKAMagicShift);
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

   delete ip_MKAMagicNumber;
   delete ip_MKAMagicShift;

   delete ip_GCWKAMagicNumber;
   delete ip_GCWKAMagicShift;

   delete ip_KATerms;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;

   delete ip_MagicKernel;
   delete ip_GCWKernel;

   xfree(il_MagicNumber);
   xfree(il_MagicShift);
   xfree(il_FactorList);
   xfree(ii_Terms);
}

void  CullenWoodallGpuWorker::CreateKernels(bool knownWorkSize)
{   
   const char *magicSource[9];
   const char *gcwSource[9];
   char  workSize[50];
   char  maxFactors[50];
   char  isCullen[50];
   char  isWoodall[50];
   char  theBase[50];

   sprintf(theBase, "#define BASE %d\n", ii_Base);
   
   if (ip_CullenWoodallApp->IsCullenSearch())
      sprintf(isCullen, "#define CHECK_CULLEN\n");
   else
      sprintf(isCullen, "\n");
   
   if (ip_CullenWoodallApp->IsWoodallSearch())
      sprintf(isWoodall, "#define CHECK_WOODALL\n");
   else
      sprintf(isWoodall, "\n");
   
   sprintf(workSize, "#define WORKSIZE %d\n", ii_WorkSize);
   sprintf(maxFactors, "#define MAX_FACTORS %d\n", ii_MaxFactorCount);
      
   magicSource[0] = workSize;
   magicSource[1] = magic;
   magicSource[2] = 0;

   gcwSource[0] = workSize;
   gcwSource[1] = theBase;
   gcwSource[2] = isCullen;
   gcwSource[3] = isWoodall;
   gcwSource[4] = maxFactors;
   gcwSource[5] = cw_kernel;
   gcwSource[6] = 0;

   ip_MagicKernel = new Kernel(ip_CullenWoodallApp->GetDevice(), "magic_kernel", magicSource);

   ip_GCWKernel = new Kernel(ip_CullenWoodallApp->GetDevice(), "cw_kernel", gcwSource);
}

void  CullenWoodallGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t group, idx;
   int32_t  n, c, ii;
   uint64_t prime;
   time_t   reportTime = time(NULL) + 60;
   
   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      ii_Groups = ip_CullenWoodallApp->GetTerms(ii_Terms, ii_MinGroupSize, ii_MaxGroupSize);
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
   }
   
   // Compute magic numbers
   ip_MagicKernel->Execute(ii_WorkSize);
      
   for (group=0; group<ii_Groups; group++)
   {
      if (ip_CullenWoodallApp->IsInterrupted() || time(NULL) > reportTime)
      {
         ip_CullenWoodallApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, group, ii_Groups);
         reportTime = time(NULL) + 60;
      }

      ii_FactorCount = 0;
      idx = group * ii_MaxGroupSize;
      ip_KATerms->SetHostMemory(&ii_Terms[idx]);

      ip_GCWKernel->Execute(ii_WorkSize);

      for (ii=0; ii<ii_FactorCount; ii++)
      {  
         idx = ii*4;
         
         n = (uint32_t) il_FactorList[idx+0];
         c = (int32_t) il_FactorList[idx+1];
         prime = il_FactorList[idx+2];
      
         if (ip_CullenWoodallApp->ReportFactor(prime, n, c))
            VerifyFactor(prime, n, c);
      }

      if (ii_FactorCount >= ii_MaxFactorCount)
         FatalError("Could not handle all GPU factors.  Sieve more deeply with the CPU before sieving with the GPU");
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
