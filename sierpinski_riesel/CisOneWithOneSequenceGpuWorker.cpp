/* CisOneWithOneSequenceGpuWorker.cpp -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CisOneWithOneSequenceGpuWorker.h"
#include "cisonesingle_kernel.gpu.h"

#define DEFAULT_HASH_MAX_DENSITY 0.65
#define HASH_MINIMUM_ELTS  8
#define HASH_MINIMUM_SHIFT 10

CisOneWithOneSequenceGpuWorker::CisOneWithOneSequenceGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ib_GpuWorker = true;
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ip_CisOneHelper = (CisOneWithOneSequenceHelper *) appHelper;
   
   ii_MaxGpuFactors = ip_SierpinskiRieselApp->GetMaxGpuFactors();
   ii_ChunksPerGpuWorker = ip_SierpinskiRieselApp->GetChunksPerGpuWorker();
}

void  CisOneWithOneSequenceGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{
   char        defines[20][50];
   const char *preKernelSources[20];
   uint32_t    defineCount = 0, idx;
   
   ii_BestQ = bestQ;
   legendre_t *legendrePtr = &ip_CisOneHelper->GetLegendre()[0];
   
   uint32_t sieveLow = ii_MinN / ii_BestQ;
   uint32_t elements = ip_CisOneHelper->GetMaxBabySteps();
   uint32_t hsize;
   
   if (HASH_MINIMUM_ELTS > elements)
      elements = HASH_MINIMUM_ELTS;
      
   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/DEFAULT_HASH_MAX_DENSITY; )
      hsize <<= 1;

   sprintf(defines[defineCount++], "#define BASE %u\n", ii_Base);
   sprintf(defines[defineCount++], "#define BESTQ %u\n", ii_BestQ);
   sprintf(defines[defineCount++], "#define SIEVE_LOW %u\n", sieveLow);
   sprintf(defines[defineCount++], "#define SEQ_PARITY %u\n", ip_FirstSequence->nParity);
   sprintf(defines[defineCount++], "#define SEQ_K %llu\n", ip_FirstSequence->k);
   sprintf(defines[defineCount++], "#define SEQ_C %lld\n", ip_FirstSequence->c);
   sprintf(defines[defineCount++], "#define KC_CORE %lld\n", ip_FirstSequence->kcCore);
   sprintf(defines[defineCount++], "#define SUBSEQUENCE_COUNT %u\n", ii_SubsequenceCount);
   sprintf(defines[defineCount++], "#define HASH_ELEMENTS %u\n", elements);
   sprintf(defines[defineCount++], "#define HASH_SIZE %u\n", hsize);
   sprintf(defines[defineCount++], "#define MAX_FACTORS %u\n", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define POWER_RESIDUE_LCM %u\n", ip_CisOneHelper->GetPowerResidueLcm());
   sprintf(defines[defineCount++], "#define DIM2 %u\n", ip_CisOneHelper->GetDim2());
   sprintf(defines[defineCount++], "#define DIM3 %u\n", ip_CisOneHelper->GetDim3());

   if (legendrePtr->haveMap)
   {
      sprintf(defines[defineCount++], "#define HAVE_LEGENDRE_TABLES");
      sprintf(defines[defineCount++], "#define LEGENDRE_MOD %u\n", legendrePtr->mod);
   }
      
   if (ip_FirstSequence->nParity == SP_MIXED)
      sprintf(defines[defineCount++], "#define HAVE_MIXED_PARITY");

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = new Kernel(ip_SierpinskiRieselApp->GetDevice(), "cisonesingle_kernel", cisonesingle_kernel, preKernelSources);

   ip_SierpinskiRieselApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_SierpinskiRieselApp->GetGpuPrimesPerWorker() * ii_ChunksPerGpuWorker;

   il_PrimeList = (uint64_t *) xmalloc(sizeof(uint64_t) * ii_WorkSize);
   
   ii_KernelWorkSize = ii_WorkSize / ii_ChunksPerGpuWorker;
   
   il_Primes = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_KernelWorkSize);
   
   if (legendrePtr->haveMap)
   {      
      if (ip_FirstSequence->nParity == SP_MIXED)
      {
         ii_DualParityMapM1 = (uint8_t *) ip_Kernel->AddCpuArgument("dualParityMapM1", sizeof(uint8_t), legendrePtr->mapSize, legendrePtr->dualParityMapM1);
         ii_DualParityMapP1 = (uint8_t *) ip_Kernel->AddCpuArgument("dualParityMapP1", sizeof(uint8_t), legendrePtr->mapSize, legendrePtr->dualParityMapP1);
      }
      else
      {
         ii_SingleParityMap = (uint8_t *) ip_Kernel->AddCpuArgument("oneParityMap", sizeof(uint8_t), legendrePtr->mapSize, legendrePtr->oneParityMap);
      }
   }
     
   ii_BabySteps   =  (uint32_t *) ip_Kernel->AddCpuArgument("babySteps", sizeof(uint32_t), ii_SubsequenceCount);
   ii_GiantSteps  =  (uint32_t *) ip_Kernel->AddCpuArgument("giantSteps", sizeof(uint32_t), ii_SubsequenceCount);

   ii_DivisorShifts       = ( int16_t *) ip_Kernel->AddCpuArgument("divisorShifts", sizeof( int16_t), ip_CisOneHelper->GetPowerResidueLcm() / 2, ip_CisOneHelper->GetDivisorShifts());
   ii_PowerResidueIndices = (uint16_t *) ip_Kernel->AddCpuArgument("powerResidueIndices", sizeof(uint32_t), ip_CisOneHelper->GetPowerResidueLcm() + 1, ip_CisOneHelper->GetPowerResidueIndices());

   ii_QIndices      = (uint32_t *) ip_Kernel->AddCpuArgument("divisorShifts", sizeof(uint32_t), ip_CisOneHelper->GetDim1(), ip_CisOneHelper->GetCongruentQIndices());
   ii_Qs            = (uint16_t *) ip_Kernel->AddCpuArgument("divisorShifts", sizeof(uint16_t), ip_CisOneHelper->GetUsedQEntries(), ip_CisOneHelper->GetAllQs());

   ii_LadderIndices = (uint32_t *) ip_Kernel->AddCpuArgument("divisorShifts", sizeof(uint32_t), ip_CisOneHelper->GetDim1(), ip_CisOneHelper->GetLadderIndices());
   ii_Ladders       = (uint16_t *) ip_Kernel->AddCpuArgument("divisorShifts", sizeof(uint16_t), ip_CisOneHelper->GetUsedLadderEntries(), ip_CisOneHelper->GetAllLadders());

   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList  = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);
   
   for (uint32_t ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      ii_BabySteps[ssIdx]  = ip_Subsequences[ssIdx].babySteps;
      ii_GiantSteps[ssIdx] = ip_Subsequences[ssIdx].giantSteps;
   }

   ip_Kernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8 + (ip_CisOneHelper->GetPowerResidueLcm()+1)*8);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CisOneWithOneSequenceGpuWorker::CleanUp(void)
{
   xfree(il_PrimeList);
   
   delete ip_Kernel;
}

void  CisOneWithOneSequenceGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx;
   uint32_t n;
   uint64_t prime;

   for (uint32_t x=0; x<ii_ChunksPerGpuWorker; x++)
   {
      memcpy(il_Primes, &il_PrimeList[ii_KernelWorkSize * x], ii_KernelWorkSize * sizeof(uint64_t));
      
      ii_FactorCount[0] = 0;

      ip_Kernel->Execute(ii_KernelWorkSize);
      
      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {
         idx = ii*2;

         n = (uint32_t) il_FactorList[idx+0];
         prime = il_FactorList[idx+1];
      
         ip_SierpinskiRieselApp->ReportFactor(prime, ip_FirstSequence, n, true);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0], ii_MaxGpuFactors);

      SetLargestPrimeTested(il_Primes[ii_KernelWorkSize-1], ii_KernelWorkSize);
   }
}

void  CisOneWithOneSequenceGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CisOneWithOneSequenceGpuWorker::TestMiniPrimeChunk not implemented");
}