/* PixSieve.cpp -- (C) Mark Rodenkirch, September 2016

   This class sets up the call to the PixSieve GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "PrimesInXGpuWorker.h"
#include "pix_kernel.gpu.h"

#define   VECTOR_SIZE   2

PrimesInXGpuWorker::PrimesInXGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    idx, defineCount;
   
   ib_GpuWorker = true;
   
   ip_PrimesInXApp = (PrimesInXApp *) theApp;
   
   ii_MaxSteps = ip_PrimesInXApp->GetMaxGpuSteps();
   ii_MaxGpuFactors = ip_PrimesInXApp->GetMaxGpuFactors() + 10;
   ii_GroupSize = (5 + ii_MaxSteps);
   
   ii_e1DigitList = BuildDigitList();
   ii_e3DigitList = BuildDigitList(1000, 3, ip_PrimesInXApp->Get3DigitTermsCopy());
   ii_e6DigitList = BuildDigitList(1000000, 6, ip_PrimesInXApp->Get6DigitTermsCopy());
   ii_e9DigitList = BuildDigitList(1000000000, 9, ip_PrimesInXApp->Get9DigitTermsCopy());
   
   defineCount = 0;
   sprintf(defines[defineCount++], "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);
   sprintf(defines[defineCount++], "#define D_EOL %u", (1 << 31));
      
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
 
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("pix_kernel", pix_kernel, preKernelSources);

   ip_PrimesInXApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_WorkSize = ip_PrimesInXApp->GetGpuPrimesPerWorker();
   
   il_PrimeList = (uint64_t *) ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t), ii_WorkSize);
   il_Residuals = (uint64_t *) ip_Kernel->AddSharedArgument("residuals", sizeof(uint64_t), ii_WorkSize);
   ii_DigitList = (uint32_t *) ip_Kernel->AddCpuArgument("digitList", sizeof(uint32_t), ii_GroupSize);
   ii_FactorCount = (uint32_t *) ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 1);
   il_FactorList = (uint64_t *) ip_Kernel->AddGpuArgument("factorList", sizeof(uint64_t), 2*ii_MaxGpuFactors);
   
   ib_Initialized = true;
}

void  PrimesInXGpuWorker::CleanUp(void)
{
   delete ip_Kernel;

   xfree(ii_e1DigitList);
   
   if (ii_e3DigitList != NULL)
      xfree(ii_e3DigitList);
   
   if (ii_e6DigitList != NULL)
      xfree(ii_e6DigitList);
   
   if (ii_e9DigitList != NULL)
      xfree(ii_e9DigitList);
}

void  PrimesInXGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t   dlIdx, termLength;
   uint64_t   prime;
   char       sPrime[20];
   uint32_t  *digitList = ii_e1DigitList;
   
   if (il_PrimeList[0] > 1000)
      digitList = ii_e3DigitList;

   if (il_PrimeList[0] > 1000000)
      digitList = ii_e6DigitList;

   if (il_PrimeList[0] > 1000000000)
      digitList = ii_e9DigitList;
   
   dlIdx = 0;
   while (digitList[dlIdx] > 0)
   {
      ii_FactorCount[0] = 0;
            
      // The first entry is the number of digits for which the residuals have been calculated.
      // The second entry is a multiplier for each term in the list.
      memcpy(ii_DigitList, &digitList[dlIdx], ii_GroupSize*sizeof(uint32_t));

      ip_Kernel->Execute(ii_WorkSize);

      for (uint32_t ii=0; ii<ii_FactorCount[0]; ii++)
      {
         uint32_t idx = ii*2;

         termLength = il_FactorList[idx+0];
         prime = il_FactorList[idx+1];
         
         sprintf(sPrime, "%llu", prime);
         
         if (strlen(sPrime) == termLength)
            ip_PrimesInXApp->ReportPrime(prime, termLength);
         else
            ip_PrimesInXApp->ReportFactor(prime, termLength);
         
         if ((ii+1) == ii_MaxGpuFactors)
            break;
      }
      
      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
      
      dlIdx += ii_GroupSize;
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  PrimesInXGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("PrimesInXGpuWorker::TestMiniPrimeChunk not implemented");
}


uint32_t  *PrimesInXGpuWorker::BuildDigitList(void)
{
   uint32_t smallTermsInList = ip_PrimesInXApp->GetMaxLength();
   uint32_t smallGroupsNeeded = 2 + (smallTermsInList / ii_MaxSteps);

   uint32_t *groupedDigitList = (uint32_t *) xmalloc(smallGroupsNeeded * ii_GroupSize * sizeof(uint32_t));
   uint32_t groupIdx = 0;
   uint32_t dlIdx = 0;
   uint32_t termLength = 0;

   uint32_t *digitList = ip_PrimesInXApp->Get1DigitTerms();
   
   for (uint32_t idx=0; idx<smallGroupsNeeded; idx++)
   {
      uint32_t groupEntries = (ii_MaxSteps < smallTermsInList ? ii_MaxSteps : smallTermsInList);
     
      groupedDigitList[groupIdx + 0] = 10;
      groupedDigitList[groupIdx + 1] = 1;
      groupedDigitList[groupIdx + 2] = termLength;
      
      memcpy(&groupedDigitList[groupIdx + 3], &digitList[dlIdx], groupEntries * sizeof(uint32_t));
      groupedDigitList[groupIdx + 3 + groupEntries] = (1 << 31);
      
      groupIdx += ii_GroupSize;
      dlIdx += groupEntries;
      termLength += groupEntries;
      
      smallTermsInList -= groupEntries;
      
      if (smallTermsInList == 0)
         break;
   }
      
   return groupedDigitList;
}

uint32_t  *PrimesInXGpuWorker::BuildDigitList(uint32_t multiplier, uint32_t power, uint32_t *digitList)
{
   uint32_t minLength = ip_PrimesInXApp->GetMinLength();
   uint32_t bigTermsInList = (minLength / power);
   uint32_t bigGroupsNeeded = bigTermsInList / ii_MaxSteps;
   
   // This ensures that all groups with big terms are "full".
   bigTermsInList = bigGroupsNeeded * ii_MaxSteps;
   
   uint32_t smallTermsInList = ip_PrimesInXApp->GetMaxLength() - (bigTermsInList * power);
   uint32_t smallGroupsNeeded = 2 + (smallTermsInList / ii_MaxSteps);

   // If I did my math correctly, this is what I expectL
   //    With these inputs:  multiplier = 1000, power = 3, minLength = 316, maxLength = 4970, maxSteps = 800
   //        we should get:  bigTermsInList = 316 /3 --> 105
   //                        bigGroupsNeeded = 105 / 800 --> 0
   //                        smallTermsInList = 4970 - (0 * 800 * 3) --> 4970
   //                        smallGroupsNeeded = (2 * 4970) / 800 --> 8
   //    This means 7 calls to the kernel to evaluate all terms as the 8th entry is empty.

   uint32_t *groupedDigitList = (uint32_t *) xmalloc((smallGroupsNeeded + bigGroupsNeeded) * ii_GroupSize * sizeof(uint32_t));
   uint32_t groupIdx = 0;
   uint32_t dlIdx = 0;
   uint32_t termLength = 0;
   
   for (uint32_t idx=0; idx<bigGroupsNeeded; idx++)
   {
      uint32_t groupEntries = (ii_MaxSteps < bigTermsInList ? ii_MaxSteps : bigTermsInList);
      
      groupedDigitList[groupIdx + 0] = multiplier;
      groupedDigitList[groupIdx + 1] = power;
      groupedDigitList[groupIdx + 2] = termLength;
      
      memcpy(&groupedDigitList[groupIdx + 3], &digitList[dlIdx], groupEntries * sizeof(uint32_t));
      groupedDigitList[groupIdx + 3 + groupEntries] = (1 << 31);
      
      groupIdx += ii_GroupSize;
      dlIdx += groupEntries;
      
      bigTermsInList -= groupEntries;
      termLength += (groupEntries * power);
      
      if (bigTermsInList == 0)
         break;
   }
   
   xfree(digitList);

   digitList = ip_PrimesInXApp->Get1DigitTerms();
   dlIdx *= power;
   
   for (uint32_t idx=0; idx<smallGroupsNeeded; idx++)
   {
      uint32_t groupEntries = (ii_MaxSteps < smallTermsInList ? ii_MaxSteps : smallTermsInList);
     
      groupedDigitList[groupIdx + 0] = 10;
      groupedDigitList[groupIdx + 1] = 1;
      groupedDigitList[groupIdx + 2] = termLength;
      
      memcpy(&groupedDigitList[groupIdx + 3], &digitList[dlIdx], groupEntries * sizeof(uint32_t));
      groupedDigitList[groupIdx + 3 + groupEntries] = (1 << 31);
      
      groupIdx += ii_GroupSize;
      dlIdx += groupEntries;
      termLength += groupEntries;
      
      smallTermsInList -= groupEntries;
      
      if (smallTermsInList == 0)
         break;
   }
  
   return groupedDigitList;
}
