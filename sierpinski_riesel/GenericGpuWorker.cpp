/* GenericGpuWorker.cpp -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "GenericGpuWorker.h"
#include "sr_kernel.h"

#define DEFAULT_HASH_MAX_DENSITY 0.65
#define HASH_MINIMUM_ELTS  8
#define HASH_MINIMUM_SHIFT 10

GenericGpuWorker::GenericGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ib_GpuWorker = true;
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);

   ii_MaxGpuFactors = ip_SierpinskiRieselApp->GetMaxGpuFactors();
}

void  GenericGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
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
   char      define10[100];
   char      define11[50];
   char      define12[50];
   char      define13[50];
   seq_t    *seq;
   uint32_t  seqIdx;
   
   ii_BestQ = bestQ;
   
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
   sprintf(define07, "#define SIEVE_RANGE %u\n", sieveRange);
   sprintf(define08, "#define SEQUENCES %u\n", ii_SequenceCount);
   sprintf(define09, "#define SUBSEQUENCES %u\n", ii_SubsequenceCount);
   sprintf(define10, "#define N_TERM(ssIdx, i, j)  ((SIEVE_LOW + (i)*BABY_STEPS + (j))*BESTQ + SUBSEQ_Q[(ssIdx)])\n"); 
   sprintf(define11, "#define HASH_ELEMENTS %u\n", elements);
   sprintf(define12, "#define HASH_SIZE %u\n", hsize);
   sprintf(define13, "#define MAX_FACTORS %u\n", ii_MaxGpuFactors);

   srSource[00] = define01;
   srSource[01] = define02;
   srSource[02] = define03;
   srSource[03] = define04;
   srSource[04] = define05;
   srSource[05] = define06;
   srSource[06] = define07;
   srSource[07] = define08;
   srSource[0x08] = define09;
   srSource[0x09] = define10;
   srSource[10] = define11;
   srSource[11] = define12;
   srSource[12] = define13;
   srSource[13] = sr_kernel;
   srSource[14] = 0;

   ip_SRKernel = new Kernel(ip_SierpinskiRieselApp->GetDevice(), "sr_kernel", srSource);
   
   AllocatePrimeList(ip_SRKernel->GetWorkGroupSize());

   il_FactorList  =  (uint64_t *) xmalloc(4*ii_MaxGpuFactors*sizeof(uint64_t));
   il_K           =  (uint64_t *) xmalloc(ii_SequenceCount*sizeof(uint64_t));
   il_C           =  (int64_t *)  xmalloc(ii_SequenceCount*sizeof(int64_t));
   ii_SeqIdx      =  (uint32_t *) xmalloc(ii_SubsequenceCount*sizeof(uint32_t));
   ii_Q           =  (uint32_t *) xmalloc(ii_SubsequenceCount*sizeof(uint32_t));
   
   seq = ip_FirstSequence;
   seqIdx = 0;
   while (seq != NULL)
   {
      il_K[seqIdx] = seq->k;
      il_C[seqIdx] = seq->c;

      for (uint32_t ssIdx=seq->ssIdxFirst; ssIdx<=seq->ssIdxLast; ssIdx++)
      {
         ii_SeqIdx[ssIdx] = seqIdx;
         ii_Q[ssIdx] = ip_Subsequences[ssIdx].q;
      }      

      seq = (seq_t *) seq->next;
      seqIdx++;
   }
   
   ip_KAPrime          = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KASeqK           = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "k", KA_HOST_TO_GPU, il_K, ii_SequenceCount);
   ip_KASeqC           = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "c", KA_HOST_TO_GPU, il_C, ii_SequenceCount);
   ip_KASubSeqSeqIdx   = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "seqIdx", KA_HOST_TO_GPU, ii_SeqIdx, ii_SubsequenceCount);
   ip_KASubSeqQ        = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "q", KA_HOST_TO_GPU, ii_Q, ii_SubsequenceCount);
   ip_KAFactorCount    = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList     = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);

   ip_SRKernel->AddArgument(ip_KAPrime);
   ip_SRKernel->AddArgument(ip_KASeqK);
   ip_SRKernel->AddArgument(ip_KASeqC);
   ip_SRKernel->AddArgument(ip_KASubSeqSeqIdx);
   ip_SRKernel->AddArgument(ip_KASubSeqQ);
   ip_SRKernel->AddArgument(ip_KAFactorCount);
   ip_SRKernel->AddArgument(ip_KAFactorList);

   ip_SRKernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8 + ii_SubsequenceCount*8);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  GenericGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KASeqK;
   delete ip_KASeqC;
   delete ip_KASubSeqSeqIdx;
   delete ip_KASubSeqQ;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;

   delete ip_SRKernel;
   
   xfree(il_FactorList);
   xfree(il_K);
   xfree(il_C);
   xfree(ii_SeqIdx);
   xfree(ii_Q);
}

void  GenericGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx, ssIdx;
   uint32_t n;
   uint64_t prime;

   ii_FactorCount = 0;

   ip_SRKernel->Execute(ii_WorkSize);

   for (uint32_t ii=0; ii<ii_FactorCount; ii++)
   {  
      idx = ii*4;
      
      ssIdx = (uint32_t) il_FactorList[idx+0];
      n = (uint32_t) il_FactorList[idx+1];
      prime = il_FactorList[idx+2];
   
      ip_SierpinskiRieselApp->ReportFactor(prime, ip_Subsequences[ssIdx].seqPtr, n, true);
      
      if (ii_FactorCount >= ii_MaxGpuFactors)
         break;
   }

   if (ii_FactorCount >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount, ii_MaxGpuFactors);

   // If we have far fewer factors than how many can fit into the buffer, shrink down the size of the buffer.
   //if (ii_FactorCount < ii_MaxGpuFactors/2)
   //{
   //   printf("shrinking from %u to %u\n", ii_MaxGpuFactors, ii_FactorCount * 2);
   //   ii_MaxGpuFactors = ii_FactorCount * 2;
   //   
   //   KernelArgument *oldArgument = ip_KAFactorList;
   //   
   //   delete oldArgument;
   //   xfree(il_FactorList);
   //   
   //   il_FactorList  =  (uint64_t *) xmalloc(4*ii_MaxGpuFactors*sizeof(uint64_t));
   //
   //   ip_KAFactorList = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, il_FactorList, 4*ii_MaxGpuFactors);
   //   
   //   ip_SRKernel->ReplaceArgument(oldArgument, ip_KAFactorList);
   //}

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  GenericGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}