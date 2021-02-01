/* CisOneWithOneSequenceGpuWorker.cpp -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "CisOneWithOneSequenceGpuWorker.h"
#include "CisOneSequenceHelper.h"
#include "cisonesingle_kernel.h"

#define DEFAULT_HASH_MAX_DENSITY 0.65
#define HASH_MINIMUM_ELTS  8
#define HASH_MINIMUM_SHIFT 10

CisOneWithOneSequenceGpuWorker::CisOneWithOneSequenceGpuWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper) : AbstractWorker(myId, theApp, appHelper)
{
   ib_GpuWorker = true;
   ip_FirstSequence = appHelper->GetFirstSequenceAndSequenceCount(ii_SequenceCount);
   ip_Subsequences = appHelper->GetSubsequences(ii_SubsequenceCount);
   
   ip_CisOneHelper = (CisOneSequenceHelper *) appHelper;

   ii_MaxGpuFactors = ip_SierpinskiRieselApp->GetMaxGpuFactors();
}

void  CisOneWithOneSequenceGpuWorker::Prepare(uint64_t largestPrimeTested, uint32_t bestQ)
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
   char      define12[50];
   char      define13[50];
   
   ii_BestQ = bestQ;
   
   uint32_t sieveLow = ii_MinN / ii_BestQ;
   uint32_t elements = ip_CisOneHelper->GetMaxBabySteps();
   uint32_t hsize;
   
   if (HASH_MINIMUM_ELTS > elements)
      elements = HASH_MINIMUM_ELTS;
      
   for (hsize = 1<<HASH_MINIMUM_SHIFT; hsize < elements/DEFAULT_HASH_MAX_DENSITY; )
      hsize <<= 1;

   sprintf(define01, "#define BASE %u\n", ii_Base);
   sprintf(define02, "#define BESTQ %u\n", ii_BestQ);
   sprintf(define03, "#define SIEVE_LOW %u\n", sieveLow);
   sprintf(define04, "#define SEQ_PARITY %u\n", ip_FirstSequence->parity);
   sprintf(define05, "#define SEQ_K %llu\n", ip_FirstSequence->k);
   sprintf(define06, "#define SEQ_C %lld\n", ip_FirstSequence->c);
   sprintf(define07, "#define KC_CORE %lld\n", ip_FirstSequence->kcCore);
   sprintf(define08, "#define SUBSEQUENCE_COUNT %u\n", ii_SubsequenceCount);
   sprintf(define09, "#define HASH_ELEMENTS %u\n", elements);
   sprintf(define10, "#define HASH_SIZE %u\n", hsize);
   sprintf(define11, "#define MAX_FACTORS %u\n", ii_MaxGpuFactors);
   sprintf(define12, "#define POWER_RESIDUE_LCM %u\n", POWER_RESIDUE_LCM);
   sprintf(define13, "#define PRL_COUNT %u\n", ip_CisOneHelper->GetPrlCount());

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
   srSource[13] = cisonesingle_kernel;
   srSource[14] = 0;

   ip_SRKernel = new Kernel(ip_SierpinskiRieselApp->GetDevice(), "cisonesingle_kernel", srSource);
   
   AllocatePrimeList(ip_SRKernel->GetWorkGroupSize());

   ip_FactorList  =  (uint64_t *) xmalloc(2*ii_MaxGpuFactors*sizeof(uint64_t));
   ip_BabySteps   =  (uint32_t *) xmalloc(ii_SubsequenceCount*sizeof(uint32_t));
   ip_GiantSteps  =  (uint32_t *) xmalloc(ii_SubsequenceCount*sizeof(uint32_t));
   
   for (uint32_t ssIdx=0; ssIdx<ii_SubsequenceCount; ssIdx++)
   {
      ip_BabySteps[ssIdx]  = ip_Subsequences[ssIdx].babySteps;
      ip_GiantSteps[ssIdx] = ip_Subsequences[ssIdx].giantSteps;
   }
   
   ip_KAPrime          = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KASubSeqBS       = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "baby_steps", KA_HOST_TO_GPU, ip_BabySteps, ii_SubsequenceCount);
   ip_KASubSeqGS       = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "giant_steps", KA_HOST_TO_GPU, ip_GiantSteps, ii_SubsequenceCount);
   ip_KADivisorShifts  = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "divisor_shifts", KA_HOST_TO_GPU, ip_CisOneHelper->GetDivisorShifts(), POWER_RESIDUE_LCM / 2);
   ip_KAPrlIndices     = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "prl_indices", KA_HOST_TO_GPU, ip_CisOneHelper->GetPrlIndices(), ip_CisOneHelper->GetPrlCount());
   
   ip_KAQIndices       = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "q_indices", KA_HOST_TO_GPU, ip_FirstSequence->congruentQIndices, ip_FirstSequence->congruentIndexCount);
   ip_KAQs             = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "qs", KA_HOST_TO_GPU, ip_FirstSequence->congruentQs, ip_FirstSequence->congruentQCount);
   ip_KALadderIndices  = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "ladder_indices", KA_HOST_TO_GPU, ip_FirstSequence->congruentLadderIndices, ip_FirstSequence->congruentIndexCount);
   ip_KALadders        = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "ladders", KA_HOST_TO_GPU, ip_FirstSequence->congruentLadders, ip_FirstSequence->congruentLadderCount);
   
   ip_KAFactorCount    = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_count", KA_BIDIRECTIONAL, &ii_FactorCount, 1);
   ip_KAFactorList     = new KernelArgument(ip_SierpinskiRieselApp->GetDevice(), "factor_list", KA_GPU_TO_HOST, ip_FactorList, 2*ii_MaxGpuFactors);

   ip_SRKernel->AddArgument(ip_KAPrime);
   ip_SRKernel->AddArgument(ip_KASubSeqBS);
   ip_SRKernel->AddArgument(ip_KASubSeqGS);
   ip_SRKernel->AddArgument(ip_KADivisorShifts);
   ip_SRKernel->AddArgument(ip_KAPrlIndices);
   ip_SRKernel->AddArgument(ip_KAQIndices);
   ip_SRKernel->AddArgument(ip_KAQs);
   ip_SRKernel->AddArgument(ip_KALadderIndices);
   ip_SRKernel->AddArgument(ip_KALadders);
   ip_SRKernel->AddArgument(ip_KAFactorCount);
   ip_SRKernel->AddArgument(ip_KAFactorList);

   ip_SRKernel->PrintStatistics(hsize * 2 + elements * 2 + (elements+1)*8 + (POWER_RESIDUE_LCM+1)*8);

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  CisOneWithOneSequenceGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KASubSeqBS;
   delete ip_KASubSeqGS;
   delete ip_KADivisorShifts;
   delete ip_KAPrlIndices;
   delete ip_KAQIndices;
   delete ip_KAQs;
   delete ip_KALadderIndices;
   delete ip_KALadders;
   delete ip_KAFactorCount;
   delete ip_KAFactorList;

   delete ip_SRKernel;
   
   xfree(ip_FactorList);
   xfree(ip_BabySteps);
   xfree(ip_GiantSteps);
}

void  CisOneWithOneSequenceGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t idx;
   uint32_t n;
   uint64_t prime;

   ii_FactorCount = 0;

   ip_SRKernel->Execute(ii_WorkSize);
   
   for (uint32_t ii=0; ii<ii_FactorCount; ii++)
   {
      idx = ii*2;

      n = (uint32_t) ip_FactorList[idx+0];
      prime = ip_FactorList[idx+1];
   
      ip_SierpinskiRieselApp->ReportFactor(prime, ip_FirstSequence, n, true);
      
      if (ii_FactorCount >= ii_MaxGpuFactors)
         break;
   }

   if (ii_FactorCount >= ii_MaxGpuFactors)
      FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factor density", ii_FactorCount, ii_MaxGpuFactors);

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  CisOneWithOneSequenceGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("CullenWoodallGpuWorker::TestMiniPrimeChunk not implemented");
}