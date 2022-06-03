/* PrimorialGpuWorker.cpp -- (C) Mark Rodenkirch, September 2016

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "PrimorialGpuWorker.h"
#include "primorial_kernel.gpu.h"

PrimorialGpuWorker::PrimorialGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;

   ib_GpuWorker = true;
   
   ip_PrimorialApp = (PrimorialApp *) theApp;

   ii_MinPrimorial = ip_PrimorialApp->GetMinPrimorial();
   ii_MaxPrimorial = ip_PrimorialApp->GetMaxPrimorial();
   
   ii_MaxGpuSteps = ip_PrimorialApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_PrimorialApp->GetMaxGpuFactors();
   
   BuildPrimeGapGroups();

   sprintf(defines[defineCount++], "#define D_IDX_OF_MIN_PRIMORIAL %d\n", ii_IdxOfMinPrimorial);
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %d\n", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define D_BIGGEST_GAP %u\n", ii_BiggestGap);
   sprintf(defines[defineCount++], "#define D_MIN_PRIMORIAL %u\n", ii_MinPrimorial);
   sprintf(defines[defineCount++], "#define FIRST_PRIMORIAL %u\n", FIRST_PRIMORIAL);
   sprintf(defines[defineCount++], "#define FIRST_PRIMORIAL_PRIME %u\n", FIRST_PRIMORIAL_PRIME);
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("primorial_kernel", primorial_kernel, preKernelSources);

   ip_PrimorialApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_PrimorialApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   il_Residues = (uint64_t *) ip_Kernel->AddSharedArgument("residues", sizeof(uint64_t), 2*ii_WorkSize);
   ii_Parameters = (uint32_t *) ip_Kernel->AddCpuArgument("parameters", sizeof(uint32_t), 2);
   ii_PrimeGaps = (uint16_t *) ip_Kernel->AddCpuArgument("primeGaps", sizeof(uint16_t), ii_MaxGpuSteps+1);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (int64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(int64_t), 4*ii_MaxGpuFactors);
   
   ip_Kernel->PrintStatistics(0);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  PrimorialGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
   
   for (uint32_t gIdx=0; gIdx<ii_GroupCount; gIdx++)
      xfree(ip_GapGroups[gIdx]);
   
   xfree(ip_GapGroups);
}

void  PrimorialGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii, idx, primorial;
   int32_t  c;
   uint64_t theFactor;

   for (uint32_t gIdx=0; gIdx<ii_GroupCount; gIdx++)
   {
      // Exit if there are no terms in this group
      if (ip_GapGroups[gIdx][0] == 0)
         break;
      
      // Set the starting index and primorial for this iteration
      ii_Parameters[0] = gIdx * ii_MaxGpuSteps;

      // Make sure we copy the 0 that marks the end of the list
      memcpy(ii_PrimeGaps, ip_GapGroups[gIdx], (1 + ii_MaxGpuSteps) * sizeof(uint16_t));
      
      ii_FactorCount[0] = 0;

      ip_Kernel->Execute(ii_WorkSize);
      
      for (ii=0; ii<ii_FactorCount[0]; ii++)
      {  
         idx = ii*4;
         
         primorial = (uint32_t) il_FactorList[idx+0];
         c = (int32_t) il_FactorList[idx+1];
         theFactor = il_FactorList[idx+2];
      
         ip_PrimorialApp->ReportFactor(theFactor, primorial, c);
            
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0], ii_MaxGpuFactors);
   }

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  PrimorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("PrimorialGpuWorker::TestMiniPrimeChunk not implemented");
}

void  PrimorialGpuWorker::BuildPrimeGapGroups(void)
{
   uint32_t *primes;
   uint16_t *primeGaps;
   uint32_t  primorialPrimeCount;
   uint32_t  remainingTerms;
   uint32_t  termsInGroup;
   
   primes = ip_PrimorialApp->GetPrimorialPrimes(primorialPrimeCount);

   for (uint32_t idx=0; idx<primorialPrimeCount; idx++)
      if (ii_MinPrimorial == primes[idx])
         ii_IdxOfMinPrimorial = idx;

   primeGaps = ip_PrimorialApp->GetPrimorialPrimeGaps(ii_BiggestGap);

   ii_GroupCount = 1 + (primorialPrimeCount / ii_MaxGpuSteps);

   ip_GapGroups = (uint16_t **) xmalloc(ii_GroupCount * sizeof(uint16_t *));

   remainingTerms = primorialPrimeCount;

   for (uint32_t gIdx=0; gIdx<ii_GroupCount; gIdx++)
   {
      termsInGroup = (ii_MaxGpuSteps > remainingTerms ? remainingTerms : ii_MaxGpuSteps);
      
      ip_GapGroups[gIdx] = (uint16_t *) xmalloc((1 + ii_MaxGpuSteps) * sizeof(uint16_t));
      memcpy(ip_GapGroups[gIdx], &primeGaps[gIdx * ii_MaxGpuSteps], termsInGroup * sizeof(uint16_t));
            
      remainingTerms -= termsInGroup;      
   }
}
