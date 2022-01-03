/* GenericGpuWorker.cpp -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "GenericGpuWorker.h"
#include "generic_kernel.h"

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
}

void  GenericGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
{ 
   ii_BestQ = bestQ;
   
   ii_KernelCount = ii_SequenceCount / ii_SequencesPerKernel;
  
   if (ii_SequenceCount % ii_SequencesPerKernel > 0)
      ii_KernelCount++;
   
   uint32_t *seqsPerKernel = (uint32_t *)  xmalloc(ii_KernelCount*sizeof(uint32_t));
   uint32_t *subseqsPerKernel = (uint32_t *)  xmalloc(ii_KernelCount*sizeof(uint32_t));
   
   il_K           = (uint64_t **) xmalloc(ii_KernelCount*sizeof(uint64_t *));
   il_C           = (int64_t **)  xmalloc(ii_KernelCount*sizeof(int64_t *));
   ii_SeqIdx      = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   ii_Q           = (uint32_t **) xmalloc(ii_KernelCount*sizeof(uint32_t *));
   ii_SubseqIdx   = (uint32_t *)  xmalloc(ii_KernelCount*sizeof(uint32_t));
   
   ip_SRKernel         = (Kernel **)   xmalloc(ii_KernelCount*sizeof(Kernel *));
   ip_KASeqK           = (KernelArgument **) xmalloc(ii_KernelCount*sizeof(KernelArgument *));
   ip_KASeqC           = (KernelArgument **) xmalloc(ii_KernelCount*sizeof(KernelArgument *));
   ip_KASubSeqSeqIdx   = (KernelArgument **) xmalloc(ii_KernelCount*sizeof(KernelArgument *));
   ip_KASubSeqQ        = (KernelArgument **) xmalloc(ii_KernelCount*sizeof(KernelArgument *));
   
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

         il_K[kIdx]       = (uint64_t *) xmalloc(seqsPerKernel[kIdx]*sizeof(uint64_t));
         il_C[kIdx]       = (int64_t *)  xmalloc(seqsPerKernel[kIdx]*sizeof(int64_t));
         ii_SeqIdx[kIdx]  = (uint32_t *) xmalloc(subseqsPerKernel[kIdx]*sizeof(uint32_t));
         ii_Q[kIdx]       = (uint32_t *) xmalloc(subseqsPerKernel[kIdx]*sizeof(uint32_t));
         
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
   
   il_FactorList  =  (uint64_t *) xmalloc(4*ii_MaxGpuFactors*sizeof(uint64_t));
   
   for (kIdx=0; kIdx<ii_KernelCount; kIdx++)
      ip_SRKernel[kIdx] = CreateKernel(kIdx, seqsPerKernel[kIdx], subseqsPerKernel[kIdx]);
   
   xfree(seqsPerKernel);
   xfree(subseqsPerKernel);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

Kernel *GenericGpuWorker::CreateKernel(uint32_t kernelIdx, uint32_t sequences, uint32_t subsequences)
{
   const char *srSource[20];
   char      define01[50];
   char      define02[50];
   char      define03[50];
   char      define04[50];
   char      define05[50];
   char      define06[50];
   char      define07[50];
   char      define08[50];
   char      define09[50];
   char      define10[50];
   char      define11[50];
   uint32_t  dIdx = 0;
  
   uint32_t r = ii_MaxN/ii_BestQ - ii_MinN/ii_BestQ + 1;
   double babyStepFactor = 1.0; // DEFAULT_BABY_STEP_FACTOR from srsieve

   uint32_t giantSteps = MAX(1, sqrt((double) r/ii_SubsequenceCount/babyStepFactor));
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

   sprintf(define01, "#define BASE %u\n", ii_Base);
   sprintf(define02, "#define BESTQ %u\n", ii_BestQ);
   sprintf(define03, "#define SIEVE_LOW %u\n", sieveLow);
   sprintf(define04, "#define SIEVE_RANGE %u\n", sieveRange);
   sprintf(define05, "#define BABY_STEPS %u\n", babySteps);
   sprintf(define06, "#define GIANT_STEPS %u\n", giantSteps);
   sprintf(define07, "#define SEQUENCES %u\n", sequences);
   sprintf(define08, "#define SUBSEQUENCES %u\n", subsequences);
   sprintf(define09, "#define HASH_ELEMENTS %u\n", elements);
   sprintf(define10, "#define HASH_SIZE %u\n", hsize);
   sprintf(define11, "#define MAX_FACTORS %u\n", ii_MaxGpuFactors);

   srSource[dIdx] = define01;
   srSource[++dIdx] = define02;
   srSource[++dIdx] = define03;
   srSource[++dIdx] = define04;
   srSource[++dIdx] = define05;
   srSource[++dIdx] = define06;
   srSource[++dIdx] = define07;
   srSource[++dIdx] = define08;
   srSource[++dIdx] = define09;
   srSource[++dIdx] = define10;
   srSource[++dIdx] = define11;
   srSource[++dIdx] = generic_kernel;
   srSource[++dIdx] = 0;

   Kernel *kernel = new Kernel(ip_SierpinskiRieselApp->GetDevice(), "generic_kernel", srSource);

   if (kernelIdx == 0)
   {
      AllocatePrimeList(kernel->GetWorkGroupSize());
      
      ip_KAPrime          = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
      ip_KAFactorCount    = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
      ip_KAFactorList     = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);
   }

   ip_KASeqK[kernelIdx]           = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "k", KA_HOST_TO_GPU, il_K[kernelIdx], sequences);
   ip_KASeqC[kernelIdx]           = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "c", KA_HOST_TO_GPU, il_C[kernelIdx], sequences);
   ip_KASubSeqSeqIdx[kernelIdx]   = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "seqIdx", KA_HOST_TO_GPU, ii_SeqIdx[kernelIdx], subsequences);
   ip_KASubSeqQ[kernelIdx]        = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "q", KA_HOST_TO_GPU, ii_Q[kernelIdx], subsequences);

   kernel->AddArgument(ip_KAPrime);
   kernel->AddArgument(ip_KASeqK[kernelIdx]);
   kernel->AddArgument(ip_KASeqC[kernelIdx]);
   kernel->AddArgument(ip_KASubSeqSeqIdx[kernelIdx]);
   kernel->AddArgument(ip_KASubSeqQ[kernelIdx]);
   kernel->AddArgument(ip_KAFactorCount);
   kernel->AddArgument(ip_KAFactorList);

   if (kernelIdx == 0)
      kernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8 + ii_SubsequenceCount*8);
   
   return kernel;
}

void  GenericGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;
   
   for (uint32_t kIdx=0; kIdx<ii_KernelCount; kIdx++)
   {
      delete ip_KASeqK[kIdx];
      delete ip_KASeqC[kIdx];
      delete ip_KASubSeqSeqIdx[kIdx];
      delete ip_KASubSeqQ[kIdx];

      delete ip_SRKernel[kIdx];

      xfree(il_K[kIdx]);
      xfree(il_C[kIdx]);
      xfree(ii_SeqIdx[kIdx]);
      xfree(ii_Q[kIdx]);
   }
   
   xfree(ip_SRKernel);
   xfree(il_FactorList);
   
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
      ii_FactorCount = 0;

      ip_SRKernel[kIdx]->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount; ii++)
      {  
         idx = ii*4;
         
         ssIdx = (uint32_t) il_FactorList[idx+0] + ii_SubseqIdx[kIdx];
         n = (uint32_t) il_FactorList[idx+1];
         prime = il_FactorList[idx+2];
      
         ip_SierpinskiRieselApp->ReportFactor(prime, ip_Subsequences[ssIdx].seqPtr, n, true);
         
         if (ii >= ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount, ii_MaxGpuFactors);
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