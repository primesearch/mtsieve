/* CisOneWithMultipleSequencesGpuWorker.cpp -- (C) Mark Rodenkirch, July 2023

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "../core/SmallHashTable.h"

#include "CisOneWithMultipleSequencesGpuWorker.h"
#include "cisonemultiple_kernel.gpu.h"

// These are from SmallHashTable.cpp
#define HASH_MINIMUM_SHIFT    8
#define HASH_MAX_DENSITY      0.60
#define HASH_MINIMUM_ELTS     8

CisOneWithMultipleSequencesGpuWorker::CisOneWithMultipleSequencesGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ib_GpuWorker = true;
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ip_CisOneHelper = (CisOneWithMultipleSequencesHelper *) appHelper;
   
   ii_MaxGpuFactors = ip_SierpinskiRieselApp->GetMaxGpuFactors();
}

void  CisOneWithMultipleSequencesGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{ 
   ii_BestQ = bestQ;

      CreateKernel();
   
   PopulateKernelArguments();
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void     CisOneWithMultipleSequencesGpuWorker::PopulateKernelArguments(void)
{
   seq_t   *seqPtr = ip_FirstSequence;
   uint32_t seqidx = 0, ssIdx = 0, seqCount = 0, subseqCount = 0;
   
   while (seqPtr != NULL)
   {
      il_Ks[seqidx]          = seqPtr->k;
      il_KCCores[seqidx]     = seqPtr->kcCore;
      
      ii_SeqData[4*seqidx+0] = seqPtr->seqIdx;
      ii_SeqData[4*seqidx+1] = seqPtr->nParity;
      ii_SeqData[4*seqidx+2] = ssIdx;
      
      for (uint32_t idx=seqPtr->ssIdxFirst; idx<=seqPtr->ssIdxLast; idx++)
      {
         ii_SubseqData[4*ssIdx+0] = ip_Subsequences[idx].q;
         ii_SubseqData[4*ssIdx+1] = ip_Subsequences[idx].babySteps;
         ii_SubseqData[4*ssIdx+2] = ip_Subsequences[idx].giantSteps;
        
         ssIdx++;
         subseqCount++;
      }
      
      ii_SeqData[4*seqidx+3] = ssIdx;
      
      seqCount++;
      seqidx++;
         
      seqPtr = (seq_t *) seqPtr->next;
   }
}

void     CisOneWithMultipleSequencesGpuWorker::CreateKernel(void)
{
   uint32_t sieveLow = ii_MinN / ii_BestQ;
   uint32_t elements = ip_CisOneHelper->GetMaxBabySteps();
   uint32_t hsize;
   
   if (HASH_MINIMUM_ELTS > elements)
      elements = HASH_MINIMUM_ELTS;
      
   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/HASH_MAX_DENSITY; )
      hsize <<= 1;

   char        defines[20][50];
   const char *preKernelSources[20];
   uint32_t    defineCount = 0, idx;

   snprintf(defines[defineCount++], 50, "#define BASE  %u", ii_Base);
   snprintf(defines[defineCount++], 50, "#define BEST_Q %u\n", ii_BestQ);
   snprintf(defines[defineCount++], 50, "#define SIEVE_LOW   %u\n", sieveLow);
   snprintf(defines[defineCount++], 50, "#define MAX_FACTORS %u\n", ii_MaxGpuFactors);
   snprintf(defines[defineCount++], 50, "#define SEQUENCES %u", ii_SequenceCount);
   snprintf(defines[defineCount++], 50, "#define SUBSEQUENCES %u", ii_SubsequenceCount);
   snprintf(defines[defineCount++], 50, "#define HASH_ELEMENTS %u\n", elements);
   snprintf(defines[defineCount++], 50, "#define HASH_SIZE %u\n", hsize);
   snprintf(defines[defineCount++], 50, "#define HASH_SIZE_M1 %u\n", hsize - 1);
   snprintf(defines[defineCount++], 50, "#define POWER_RESIDUE_LCM %u\n", ip_CisOneHelper->GetPowerResidueLcm());
   snprintf(defines[defineCount++], 50, "#define DIM2 %u\n", ip_CisOneHelper->GetDim2());
   snprintf(defines[defineCount++], 50, "#define DIM3 %u\n", ip_CisOneHelper->GetDim3());
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("cisonemultiple_kernel", cisonemultiple_kernel, preKernelSources);

   ip_SierpinskiRieselApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_SierpinskiRieselApp->GetGpuPrimesPerWorker();

   il_Primes      = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_PrimesInList);
   
   il_Ks          = (uint64_t *) ip_Kernel->AddCpuArgument("ks", sizeof(uint64_t), ii_SequenceCount);
   il_KCCores     = ( int64_t *) ip_Kernel->AddCpuArgument("kccores", sizeof(int64_t), ii_SequenceCount);
   ii_SeqData     = (uint32_t *) ip_Kernel->AddCpuArgument("seqdata", sizeof(uint32_t), 4 * (ii_SequenceCount));
   ii_SubseqData  = (uint32_t *) ip_Kernel->AddCpuArgument("subseqdata", sizeof(uint32_t), 4 * ii_SubsequenceCount);
   
   ii_DivisorShifts          = ( int16_t *) ip_Kernel->AddCpuArgument("divisorShifts", sizeof( int16_t), ip_CisOneHelper->GetPowerResidueLcm() / 2, ip_CisOneHelper->GetDivisorShifts());
   ii_PowerResidueIndices    = (uint16_t *) ip_Kernel->AddCpuArgument("powerResidueIndices", sizeof(uint16_t), ip_CisOneHelper->GetPowerResidueLcm() + 1, ip_CisOneHelper->GetPowerResidueIndices());
   ii_CongruentSubseqIndices = (uint32_t *) ip_Kernel->AddCpuArgument("congruentSubseqIndices", sizeof(uint32_t), ip_CisOneHelper->GetDim1(), ip_CisOneHelper->GetCongruentSubseqIndices());
   ii_CongruentSubseqs       = (uint32_t *) ip_Kernel->AddCpuArgument("conguentSubseqs", sizeof(uint32_t), ip_CisOneHelper->GetUsedSubseqEntries(), ip_CisOneHelper->GetAllSubseqs());

   ii_Ladders     = (uint16_t *) ip_Kernel->AddCpuArgument("ladders", sizeof(uint16_t), ip_CisOneHelper->GetUsedLadderEntries(), ip_CisOneHelper->GetAllLadders());
   
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList  = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8);
}

void  CisOneWithMultipleSequencesGpuWorker::CleanUp(void)
{
   xfree(ip_Kernel);
   
   xfree(il_Primes);
   xfree(il_Ks);
   xfree(il_KCCores);
   xfree(ii_SeqData);
   xfree(ii_SubseqData);
   xfree(ii_DivisorShifts);
   xfree(ii_PowerResidueIndices);
   xfree(ii_CongruentSubseqIndices);
   xfree(ii_CongruentSubseqs);
   xfree(ii_Ladders);  
   xfree(ii_FactorCount);
   xfree(il_FactorList);
}

void  CisOneWithMultipleSequencesGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, ssIdx;
   uint32_t n, gpuFactors;
   uint64_t prime;

   // All kernels share the same memory for the primes
   memcpy(il_Primes, il_PrimeList, sizeof(uint64_t) * ii_PrimesInList);

   ii_FactorCount[0] = 0;

   ip_Kernel->Execute(ii_PrimesInList);

   gpuFactors = ii_FactorCount[0];
   
   for (uint32_t ii=0; ii<gpuFactors; ii++)
   {  
      idx = ii*4;
      
      prime = il_FactorList[idx+0];
      ssIdx = (uint32_t) il_FactorList[idx+1];
      n = (uint32_t) il_FactorList[idx+2];

      if ((ii+1) == ii_MaxGpuFactors)
         break;
   
      //printf("%llu %llu %llu %llu\n", il_FactorList[idx+0], il_FactorList[idx+1], il_FactorList[idx+2], il_FactorList[idx+3]);
      ip_SierpinskiRieselApp->ReportFactor(prime, ip_Subsequences[ssIdx].seqPtr, n, true);
      
   }

   if (gpuFactors >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", gpuFactors, ii_MaxGpuFactors);

   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  CisOneWithMultipleSequencesGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}
