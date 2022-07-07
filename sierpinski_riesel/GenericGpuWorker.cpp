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

#define PARAM_SEQUENCE        0
#define PARAM_SUBSEQUENCE     1
#define PARAM_SIEVE_RANGE     2
#define PARAM_BABY_STEPS      3
#define PARAM_GIANT_STEPS     4

// These are only used by the Worker.  They are not used by the Kernel.
#define PARAM_SUBSEQ_IDX     11
#define PARAM_HASH_SIZE      12
#define PARAM_HASH_ELEMENTS  13

#define MAX_PARAMETERS       19

GenericGpuWorker::GenericGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ib_GpuWorker = true;
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ib_CanUseCIsOneLogic = ip_SierpinskiRieselApp->CanUseCIsOneLogic();
   il_MaxK = ip_SierpinskiRieselApp->GetMaxK();
   ib_HaveSingleC = ip_SierpinskiRieselApp->HasSingleC();

   ii_MaxGpuFactors = ip_SierpinskiRieselApp->GetMaxGpuFactors();
   ii_KernelCount = ip_SierpinskiRieselApp->GetKernelCount();
   ii_ChunksPerGpuWorker = ip_SierpinskiRieselApp->GetChunksPerGpuWorker();
   
   il_PrimeList = NULL;
}

void  GenericGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{ 
   ii_BestQ = bestQ;

   uint32_t sequencesPerKernel = ii_SequenceCount / ii_KernelCount;
   
   while (sequencesPerKernel * ii_KernelCount < ii_SequenceCount)
      sequencesPerKernel++; 
   
   CreateKernels(sequencesPerKernel);
   
   PopulateKernelArguments(sequencesPerKernel);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void    GenericGpuWorker::CreateKernels(uint32_t sequencesPerKernel)
{
   seq_t   *seqPtr = ip_FirstSequence;
   uint32_t kIdx = 0;
   uint32_t seqCount = 0;
   uint32_t ssCount = 0;
   
   ii_SubseqIdx   = (uint32_t  *) xmalloc(ii_KernelCount*sizeof(uint32_t));
   ip_Kernel      = (GpuKernel **) xmalloc(ii_KernelCount*sizeof(GpuKernel *));
   
   il_Primes      = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   il_Ks          = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   il_Cs          = ( int64_t **) xmalloc(ii_KernelCount*sizeof( int64_t *));
   ii_SeqIdxs     = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   ii_Qs          = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   ii_FactorCount = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   il_FactorList  = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   
   ii_SubseqIdx[kIdx] = seqPtr->ssIdxFirst;
   
   // Create the Kernels which in turn will allocate memory for the kernel arguments
   while (seqPtr != NULL)
   {
      if (seqCount == sequencesPerKernel)
      {
         ip_Kernel[kIdx] = CreateKernel(kIdx, seqCount, ssCount);
         
         kIdx++;
         
         seqCount = 0;
         ssCount = 0;
         ii_SubseqIdx[kIdx] = seqPtr->ssIdxFirst;
      }
      
      seqCount++;
      ssCount += seqPtr->ssCount;

      seqPtr = (seq_t *) seqPtr->next;
   }
         
   // Create the last Kernel
   ip_Kernel[kIdx] = CreateKernel(kIdx, seqCount, ssCount);
}

void    GenericGpuWorker::PopulateKernelArguments(uint32_t sequencesPerKernel)
{
   seq_t   *seqPtr = ip_FirstSequence;
   uint32_t kIdx = 0;
   uint32_t kSeqIdx = 0;
   uint32_t kSubseqIdx = 0;
      
   // Now we can populate the memory for the kernel arguments
   while (seqPtr != NULL)
   {
      if (kSeqIdx == sequencesPerKernel)
      {         
         kIdx++;
         
         kSeqIdx = 0;
         kSubseqIdx = 0;
      }
      
      il_Ks[kIdx][kSeqIdx] = seqPtr->k;
      il_Cs[kIdx][kSeqIdx] = seqPtr->c;

      for (uint32_t ssIdx=seqPtr->ssIdxFirst; ssIdx<=seqPtr->ssIdxLast; ssIdx++)
      {
         ii_SeqIdxs[kIdx][kSubseqIdx] = kSeqIdx;
         ii_Qs[kIdx][kSubseqIdx] = ip_Subsequences[ssIdx].q;
         
         kSubseqIdx++;
      }      

      kSeqIdx++;
      seqPtr = (seq_t *) seqPtr->next;
   }   
}

GpuKernel *GenericGpuWorker::CreateKernel(uint32_t kIdx, uint32_t sequences, uint32_t subsequences)
{ 
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

   char        defines[20][50];
   const char *preKernelSources[20];
   uint32_t    defineCount = 0, idx;

   sprintf(defines[defineCount++], "#define BASE  %u", ii_Base);
   sprintf(defines[defineCount++], "#define BESTQ %u", ii_BestQ);
   sprintf(defines[defineCount++], "#define SIEVE_LOW   %u", sieveLow);
   sprintf(defines[defineCount++], "#define SIEVE_RANGE %u", sieveRange);
   sprintf(defines[defineCount++], "#define MAX_FACTORS %u", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define BABY_STEPS %u", babySteps);
   sprintf(defines[defineCount++], "#define GIANT_STEPS %u", giantSteps);
   sprintf(defines[defineCount++], "#define SEQUENCES %u", sequences);
   sprintf(defines[defineCount++], "#define SUBSEQUENCES %u", subsequences);
   sprintf(defines[defineCount++], "#define HASH_ELEMENTS %u", elements);
   sprintf(defines[defineCount++], "#define HASH_SIZE %u", hsize);
   
   if (ib_HaveSingleC)
      sprintf(defines[defineCount++], "#define SINGLE_C %" PRId64"", ip_FirstSequence->c);
   
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;

   GpuKernel *kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("generic_kernel", generic_kernel, preKernelSources);

   ip_SierpinskiRieselApp->SetGpuWorkGroupSize(kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_SierpinskiRieselApp->GetGpuPrimesPerWorker() * ii_ChunksPerGpuWorker;
   ii_KernelWorkSize = ii_PrimesInList / ii_ChunksPerGpuWorker;

   if (il_PrimeList == NULL)
      il_PrimeList      = (uint64_t *) xmalloc(sizeof(uint64_t) * ii_PrimesInList);

   il_Primes[kIdx]      = (uint64_t *) kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_KernelWorkSize);
   il_Ks[kIdx]          = (uint64_t *) kernel->AddCpuArgument("k", sizeof(uint64_t), sequences);
   il_Cs[kIdx]          = ( int64_t *) kernel->AddCpuArgument("c", sizeof(int64_t), sequences);
   ii_SeqIdxs[kIdx]     = (uint32_t *) kernel->AddCpuArgument("seqIdx", sizeof(uint32_t), subsequences);
   ii_Qs[kIdx]          = (uint32_t *) kernel->AddCpuArgument("q", sizeof(uint32_t), subsequences);
   ii_FactorCount[kIdx] = (uint32_t *) kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList[kIdx]  = (uint64_t *) kernel->AddGpuArgument("factorList", sizeof(uint64_t), 4*ii_MaxGpuFactors);

   if (kIdx == 0)
      kernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8 + ii_SubsequenceCount*8);
   
   return kernel;
}

void  GenericGpuWorker::CleanUp(void)
{ 
   for (uint32_t kIdx=0; kIdx<ii_KernelCount; kIdx++)
      delete ip_Kernel[kIdx];
   
   xfree(ip_Kernel);
   xfree(il_Ks);
   xfree(il_Cs);  
   xfree(ii_SeqIdxs);  
   xfree(ii_Qs);  
   xfree(ii_FactorCount);
   xfree(il_FactorList);
   
   xfree(ii_SubseqIdx);
   xfree(il_PrimeList);
}

void  GenericGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, ssIdx;
   uint32_t n;
   uint64_t prime;

   for (uint32_t cIdx=0; cIdx<ii_ChunksPerGpuWorker; cIdx++)
   {
      for (uint32_t kIdx=0; kIdx<ii_KernelCount; kIdx++)
      {
         uint32_t pIndex = ii_KernelWorkSize * cIdx;
         
         memcpy(il_Primes[kIdx], &il_PrimeList[pIndex], sizeof(uint64_t) * ii_KernelWorkSize);
         
         ii_FactorCount[kIdx][0] = 0;

         ip_Kernel[kIdx]->Execute(ii_KernelWorkSize);

         for (uint32_t ii=0; ii<ii_FactorCount[kIdx][0]; ii++)
         {  
            idx = ii*4;
            
            ssIdx = (uint32_t) il_FactorList[kIdx][idx+0] + ii_SubseqIdx[kIdx];
            n = (uint32_t) il_FactorList[kIdx][idx+1];
            prime = il_FactorList[kIdx][idx+2];

            if ((ii+1) == ii_MaxGpuFactors)
               break;
         
            ip_SierpinskiRieselApp->ReportFactor(prime, ip_Subsequences[ssIdx].seqPtr, n, true);
            
         }

         //if (ii_FactorCount[kIdx][0] >= ii_MaxGpuFactors)
         //   FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount[kIdx][0], ii_MaxGpuFactors);
      }
   }

   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
   
   // Determine if we can switch to the CisOne workers.  This will automatically switch
   // to the CisOne GPU workers.
   if (ib_CanUseCIsOneLogic && il_PrimeList[ii_PrimesInList-1] > il_MaxK && ii_SequenceCount == 1)
      ip_SierpinskiRieselApp->SetRebuildNeeded();
}

void  GenericGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}