/* GenericGpuWorker.cpp -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "GenericGpuWorker.h"
#include "generic_kernel.gpu.h"

#define DEFAULT_HASH_MAX_DENSITY 0.65
#define HASH_MINIMUM_ELTS  8
#define HASH_MINIMUM_SHIFT 10

GenericGpuWorker::GenericGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ib_GpuWorker = true;
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ib_CanUseCIsOneLogic = ip_SierpinskiRieselApp->CanUseCIsOneLogic();
   il_MaxK = ip_SierpinskiRieselApp->GetMaxK();

   ii_MaxGpuFactors = ip_SierpinskiRieselApp->GetMaxGpuFactors();
   ii_SequencesPerKernel = ip_SierpinskiRieselApp->GetSequecesPerKernel();
   
   // We need to explicitly allocate in this GPU worker since there are multiple Kernels, each with its own list
   il_PrimeList = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
}

void  GenericGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{ 
   ii_BestQ = bestQ;
   
   ii_KernelCount = ii_SequenceCount / ii_SequencesPerKernel;
  
   if (ii_SequenceCount % ii_SequencesPerKernel > 0)
      ii_KernelCount++;
   
   uint32_t *seqsPerKernel = (uint32_t *)  xmalloc(ii_KernelCount*sizeof(uint32_t));
   uint32_t *subseqsPerKernel = (uint32_t *)  xmalloc(ii_KernelCount*sizeof(uint32_t));
   
   ip_Kernel      = (Kernel   **) xmalloc(ii_KernelCount*sizeof(Kernel));
   il_Primes      = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   il_K           = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   il_C           = ( int64_t **) xmalloc(ii_KernelCount*sizeof( int64_t *));
   ii_SeqIdx      = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   ii_Q           = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   ii_SubseqIdx   = (uint32_t  *) xmalloc(ii_KernelCount*sizeof(uint32_t  ));
   ii_FactorCount = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   il_FactorList  = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   
   seq_t   *seqPtr = ip_FirstSequence;
   uint32_t seqIdx = 0;
   uint32_t kIdx = 0;
   uint32_t kSubseqIdx = 0;
   
   while (seqPtr != NULL)
   {
      if (seqIdx == ii_SequencesPerKernel)
      {
         kIdx++;
         seqIdx = 0;
      }
      
      seqsPerKernel[kIdx]++;
      subseqsPerKernel[kIdx] += seqPtr->ssCount;

      seqPtr = (seq_t *) seqPtr->next;
      seqIdx++;
   }
   
   seqPtr = ip_FirstSequence;
   seqIdx = 0;
   kIdx = 0;
   ii_SubseqIdx[0] = 0;
   
   while (seqPtr != NULL)
   {
      if (seqPtr == ip_FirstSequence || seqIdx == ii_SequencesPerKernel)
      {
         if (seqIdx == ii_SequencesPerKernel)
            kIdx++;

         ip_Kernel[kIdx] = CreateKernel(kIdx, seqsPerKernel[kIdx], subseqsPerKernel[kIdx]);
                  
         seqIdx = 0;
         kSubseqIdx = 0;
      }
      
      il_K[kIdx][seqIdx] = seqPtr->k;
      il_C[kIdx][seqIdx] = seqPtr->c;

      if (seqIdx == 0)
         ii_SubseqIdx[kIdx] = seqPtr->ssIdxFirst;

      for (uint32_t ssIdx=seqPtr->ssIdxFirst; ssIdx<=seqPtr->ssIdxLast; ssIdx++)
      {
         ii_SeqIdx[kIdx][kSubseqIdx] = seqIdx;
         ii_Q[kIdx][kSubseqIdx] = ip_Subsequences[ssIdx].q;
         kSubseqIdx++;
      }      

      seqPtr = (seq_t *) seqPtr->next;
      seqIdx++;
   }

   xfree(seqsPerKernel);
   xfree(subseqsPerKernel);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

Kernel *GenericGpuWorker::CreateKernel(uint32_t kIdx, uint32_t sequences, uint32_t subsequences)
{
   char        defines[20][50];
   const char *preKernelSources[20];
   uint32_t    defineCount = 0, idx;
  
   uint32_t r = ii_MaxN/ii_BestQ - ii_MinN/ii_BestQ + 1;
   double babyStepFactor = 1.0; // DEFAULT_BABY_STEP_FACTOR from srsieve

   uint32_t giantSteps = MAX(1, sqrt((double) r/subsequences/babyStepFactor));
   uint32_t babySteps = MIN(r, ceil((double) r/giantSteps));

   if (babySteps > HASH_MAX_ELTS)
   {
      giantSteps = ceil((double)r/HASH_MAX_ELTS);
      babySteps = ceil((double)r/giantSteps);
   }
   
   uint32_t sieveLow = ii_MinN / ii_BestQ;
   uint32_t sieveRange = babySteps * giantSteps;
   uint32_t elements = babySteps;
   uint32_t hsize;
   
   if (HASH_MINIMUM_ELTS > elements)
      elements = HASH_MINIMUM_ELTS;
      
   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/DEFAULT_HASH_MAX_DENSITY; )
      hsize <<= 1;

   sprintf(defines[defineCount++], "#define BASE %u", ii_Base);
   sprintf(defines[defineCount++], "#define BESTQ %u", ii_BestQ);
   sprintf(defines[defineCount++], "#define SIEVE_LOW %u", sieveLow);
   sprintf(defines[defineCount++], "#define SIEVE_RANGE %u", sieveRange);
   sprintf(defines[defineCount++], "#define BABY_STEPS %u", babySteps);
   sprintf(defines[defineCount++], "#define GIANT_STEPS %u", giantSteps);
   sprintf(defines[defineCount++], "#define SEQUENCES %u", sequences);
   sprintf(defines[defineCount++], "#define SUBSEQUENCES %u", subsequences);
   sprintf(defines[defineCount++], "#define HASH_ELEMENTS %u", elements);
   sprintf(defines[defineCount++], "#define HASH_SIZE %u", hsize);
   sprintf(defines[defineCount++], "#define MAX_FACTORS %u", ii_MaxGpuFactors);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   Kernel *kernel = new Kernel(ip_SierpinskiRieselApp->GetDevice(), "generic_kernel", generic_kernel, preKernelSources);

   ip_SierpinskiRieselApp->SetGpuWorkGroupSize(kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_SierpinskiRieselApp->GetGpuPrimesPerWorker();

   il_Primes[kIdx]       = (uint64_t *) kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   il_K[kIdx]            = (uint64_t *) kernel->AddGpuArgument("k", sizeof(uint64_t), sequences);
   il_C[kIdx]            = ( int64_t *) kernel->AddCpuArgument("c", sizeof(int64_t), sequences);
   ii_SeqIdx[kIdx]       = (uint32_t *) kernel->AddCpuArgument("seqidx", sizeof(uint32_t), subsequences);
   ii_Q[kIdx]            = (uint32_t *) kernel->AddCpuArgument("q", sizeof(uint32_t), subsequences);
   ii_FactorCount[kIdx]  = (uint32_t *) kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList[kIdx]   = (uint64_t *) kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   if (kIdx == 0)
      kernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8 + ii_SubsequenceCount*8);
   
   return kernel;
}

void  GenericGpuWorker::CleanUp(void)
{
   for (uint32_t kIdx=0; kIdx<ii_KernelCount; kIdx++)
      delete ip_Kernel[kIdx];
   
   delete ip_Kernel;
   
   xfree(il_PrimeList);
   
   xfree(il_Primes);
   xfree(il_K);
   xfree(il_C);
   xfree(ii_SeqIdx);
   xfree(ii_Q);
   xfree(ii_SubseqIdx);
}

void  GenericGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, ssIdx;
   uint32_t n;
   uint64_t prime;

   for (uint32_t kIdx=0; kIdx<ii_KernelCount; kIdx++)
   {
      memcpy(il_Primes[kIdx], il_PrimeList, ii_WorkSize * sizeof(uint64_t));
      
      ii_FactorCount[0][0] = 0;

      ip_Kernel[kIdx]->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount[0][0]; ii++)
      {  
         idx = ii*4;
         
         ssIdx = (uint32_t) il_FactorList[kIdx][idx+0] + ii_SubseqIdx[kIdx];
         n = (uint32_t) il_FactorList[kIdx][idx+1];
         prime = il_FactorList[kIdx][idx+2];
      
         ip_SierpinskiRieselApp->ReportFactor(prime, ip_Subsequences[ssIdx].seqPtr, n, true);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0][0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[0][0], ii_MaxGpuFactors);
   }

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
   
   // Determine if we can switch to the CisOne workers.  This will automatically switch
   // to the CisOne GPU workers.
   if (ib_CanUseCIsOneLogic && il_PrimeList[ii_WorkSize-1] > il_MaxK && ii_SequenceCount == 1)
      ip_SierpinskiRieselApp->SetRebuildNeeded();
}

void  GenericGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}