/* CarolKyneaGpuWorker.cpp -- (C) Mark Rodenkirch, November 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include <assert.h>

#include "../core/HashTable.h"
#include "CarolKyneaGpuWorker.h"
#include "ck_kernel.gpu.h"

#define DEFAULT_HASH_MAX_DENSITY 0.65
#define HASH_MINIMUM_ELTS     8
#define HASH_MINIMUM_SHIFT   11

#define ROOT_COUNT 4

CarolKyneaGpuWorker::CarolKyneaGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[20][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   uint32_t    giantSteps, babySteps, sieveLow, sieveRange;
   
   ib_GpuWorker = true;
   ip_CarolKyneaApp = (CarolKyneaApp *) theApp;
      
   ii_Base = ip_CarolKyneaApp->GetBase();
   ii_MinN = ip_CarolKyneaApp->GetMinN();
   ii_MaxN = ip_CarolKyneaApp->GetMaxN();
   ii_MaxGpuFactors = ip_CarolKyneaApp->GetMaxGpuFactors();

   uint32_t r = ii_MaxN - ii_MinN + 1;

   // In the worst case we will do do one table insertion and one mulmod
   // for m baby steps, then s table lookups and s mulmods for M giant
   // steps. The average case depends on how many solutions are found
   // and how early in the loop they are found, which I don't know how
   // to analyse. However for the worst case we just want to minimise
   // m + s*M subject to m*M >= r, which is when m = sqrt(s*r).
  
   giantSteps = MAX(1, sqrt((double) r/ROOT_COUNT));
   babySteps = MIN(r, ceil((double) r/giantSteps));

   if (babySteps > HASH_MAX_ELTS)
   {
      giantSteps = ceil((double)r/HASH_MAX_ELTS);
      babySteps = ceil((double)r/giantSteps);
   }

   sieveLow = ii_MinN;
   sieveRange = babySteps*giantSteps;
   
   assert(sieveLow <= ip_CarolKyneaApp->GetMinN());
   assert(ii_MaxN < sieveLow+sieveRange);

   uint32_t elements = babySteps;
   uint32_t hsize;
   
   if (HASH_MINIMUM_ELTS > elements)
      elements = HASH_MINIMUM_ELTS;
      
   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/DEFAULT_HASH_MAX_DENSITY; )
      hsize <<= 1;

   sprintf(defines[defineCount++], "#define BASE    %u", ii_Base);
   sprintf(defines[defineCount++], "#define N_MIN   %u", ii_MinN);
   sprintf(defines[defineCount++], "#define SIEVE_LOW   %u", sieveLow);
   sprintf(defines[defineCount++], "#define SIEVE_RANGE %u", sieveRange);
   sprintf(defines[defineCount++], "#define MAX_FACTORS %u", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define BABY_STEPS %u", babySteps);
   sprintf(defines[defineCount++], "#define GIANT_STEPS %u", giantSteps);
   sprintf(defines[defineCount++], "#define HASH_ELEMENTS %u", elements);
   sprintf(defines[defineCount++], "#define HASH_SIZE %u", hsize);
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("ck_kernel", ck_kernel, preKernelSources);

   ip_CarolKyneaApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_CarolKyneaApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (int64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(int64_t), 4*ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(0);
      
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CarolKyneaGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  CarolKyneaGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx;
   int32_t  n, c;
   uint64_t prime;
  
   ii_FactorCount[0] = 0;
   
   ip_Kernel->Execute(ii_PrimesInList);

   for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
   {  
      idx = ii*4;
      
      n = (uint32_t) il_FactorList[idx+0];
      c = (int32_t) il_FactorList[idx+1];
      prime = il_FactorList[idx+2];

      ip_CarolKyneaApp->ReportFactor(prime, n, c);
      
      if ((ii+1) == ii_MaxGpuFactors)
         break;
   }

   if (ii_FactorCount[0] >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);

   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  CarolKyneaGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CarolKyneaGpuWorker::TestMiniPrimeChunk not implemented");
}
