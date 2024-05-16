/* LifchitzGpuWorker.cpp -- (C) Mark Rodenkirch, April 2024

   This class sets up the call to the LifchitzGpuWorker GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "LifchitzGpuWorker.h"
#include "lifchitz_kernel.gpu.h"

#define BIT(base)       (base - ii_MinBase)

LifchitzGpuWorker::LifchitzGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{ 
   char        defines[10][50];
   const char *preKernelSources[10];
   uint32_t    defineCount = 0, idx;
   
   ib_GpuWorker = true;
   
   ip_LifchitzApp = (LifchitzApp *) theApp;

   ii_MinX = ip_LifchitzApp->GetMinX();
   ii_MaxX = ip_LifchitzApp->GetMaxX();
   
   ii_MaxGpuFactors = ip_LifchitzApp->GetMaxGpuFactors();

   AllocateMemory();

   snprintf(defines[defineCount++], 50, "#define D_MAX_X_BASES %u", ii_XPerChunk);
   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);

   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("lifchitz_kernel", lifchitz_kernel, preKernelSources);

   ip_LifchitzApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_LifchitzApp->GetGpuPrimesPerWorker();
   
   il_PrimeList    = (uint64_t *)  ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t),  ii_PrimesInList);
   ii_KernelXBases = (uint32_t *)  ip_Kernel->AddCpuArgument("xbases", sizeof(uint32_t),  ii_XPerChunk);
   ip_TermList     = (gputerm_t *) ip_Kernel->AddCpuArgument("terms",  sizeof(gputerm_t), ii_MaxTermsPerGroup);
   //ip_SignList     = (uint8_t *)   ip_Kernel->AddCpuArgument("signs",  sizeof(uint8_t),   ii_MaxTermsPerGroup);
   ii_FactorCount  = (uint32_t *)  ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 10);
   ip_FactorList   = (factor_t *)  ip_Kernel->AddGpuArgument("factorList", sizeof(factor_t), ii_MaxGpuFactors);

   ip_Kernel->PrintStatistics(0);

   // The thread can't start until initialization is done
   ib_Initialized = true;

   il_NextTermsBuild = 0;
}

void  LifchitzGpuWorker::CleanUp(void)
{
   delete ip_Kernel;
}

void  LifchitzGpuWorker::AllocateMemory(void)
{
   uint32_t idx, xIdx;

   ii_XPerChunk = ip_LifchitzApp->GetXPerChunk();

   uint32_t minXForChunk = ii_MinX;
   
   ii_XChunks = 0;
   while (minXForChunk < ii_MaxX)
   {
      ii_XChunks++;
      minXForChunk += ii_XPerChunk;
   }
   
   if (ii_XPerChunk > ii_MaxX - ii_MinX)
   {
      ii_XChunks = 1;
      ii_XPerChunk = 1 + ii_MaxX - ii_MinX;
   }
   
   ip_XBases = (uint32_t **) xmalloc(ii_XPerChunk, sizeof(uint32_t *), "xBases *");

   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
      ip_XBases[xIdx] = (uint32_t *) xmalloc(ii_XPerChunk, sizeof(uint32_t), "xBases");
   
   ip_TermsInGroup = (uint32_t *) xmalloc(ii_XChunks, sizeof(uint32_t *), "termsInGroup");
      
   term_t *terms = ip_LifchitzApp->GetTerms();
   idx = 0;

   while (terms[idx].x > 0)
   {
      xIdx = (terms[idx].x - ii_MinX) / ii_XPerChunk;
      
      ip_TermsInGroup[xIdx]++;

      idx++;
   }

   if (ip_LifchitzApp->GetGpuWorkerCount() + ip_LifchitzApp->GetCpuWorkerCount() > 1)
      xfree(terms);
   
   ii_MaxTermsPerGroup = 0; 
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {       
      if (ip_TermsInGroup[xIdx] > ii_MaxTermsPerGroup)
            ii_MaxTermsPerGroup = ip_TermsInGroup[xIdx];
   }
      
   ip_GpuTerms    = (gputerm_t **) xmalloc(ii_XChunks, sizeof(gputerm_t *), "gpuTerms *");
   ip_GpuSigns    = (uint8_t **)   xmalloc(ii_XChunks, sizeof(uint8_t *),   "gpuSigns *");
   ip_TermIndexes = (uint64_t **)  xmalloc(ii_XChunks, sizeof(uint64_t *),  "termIndexes *");

   minXForChunk = ii_MinX;
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      ip_GpuTerms[xIdx]    = (gputerm_t *) xmalloc(1 + ii_MaxTermsPerGroup, sizeof(gputerm_t), "gpuTerms");
      ip_GpuSigns[xIdx]    = (uint8_t *)   xmalloc(1 + ii_MaxTermsPerGroup, sizeof(uint8_t),   "gpuSigns");
      ip_TermIndexes[xIdx] = (uint64_t *)  xmalloc(1 + ii_MaxTermsPerGroup, sizeof(uint64_t),  "termIndexes");

      minXForChunk += ii_XPerChunk;
   }
}

void  LifchitzGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t  xIdx;
   time_t    reportTime = time(NULL) + 60;

   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      CreateTermGroups();
      
      il_NextTermsBuild = (il_PrimeList[ii_PrimesInList-1] << 1);
   }
         
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      if (ip_LifchitzApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_LifchitzApp->WriteToConsole(COT_SIEVE, "Thread %d was interrupted.  Completed %u of %u iterations", ii_MyId, xIdx, ii_XChunks);
         reportTime = time(NULL) + 60;
      }
      
      memcpy(ii_KernelXBases, ip_XBases[xIdx],   ii_XPerChunk        * sizeof(uint32_t));
      memcpy(ip_TermList,     ip_GpuTerms[xIdx], ii_MaxTermsPerGroup * sizeof(gputerm_t));
     // memcpy(ip_SignList,     ip_GpuSigns[xIdx], ii_MaxTermsPerGroup * sizeof(uint8_t));

      ii_FactorCount[0] = 0;
      
      ip_Kernel->Execute(ii_PrimesInList);
      
      for (uint32_t fIdx=0; fIdx<ii_FactorCount[0]; fIdx++)
      {
         uint64_t tIdx =  ip_FactorList[fIdx].termIdx;
         uint32_t x = (uint32_t) ip_TermList[tIdx].x;
         uint32_t y = (uint32_t) ip_TermList[tIdx].y;
         int32_t  sign = (ip_FactorList[fIdx].sign == 1 ? +1 : -1);
         uint64_t theFactor = ip_FactorList[fIdx].factor;
         
         ip_LifchitzApp->ReportFactor(theFactor, x, y, sign, ip_TermIndexes[xIdx][tIdx]);
         
         if ((fIdx+1) == ii_MaxGpuFactors)
            break;
      }

      if (ii_FactorCount[0] >= ii_MaxGpuFactors)
         FatalError("Could not handle all GPU factors.  A range of p generated %u factors (limited to %u).  Use -M to increase max factors", ii_FactorCount[0], ii_MaxGpuFactors);
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_PrimesInList-1], ii_PrimesInList);
}

void  LifchitzGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("LifchitzGpuWorker::TestMiniPrimeChunk not implemented");
}

void   LifchitzGpuWorker::CreateTermGroups(void)
{
   uint32_t idx, xIdx, tIdx;
   uint32_t xOffset;

   // These arrays will have zeros where we don't need to compute x^x or y^y
   uint32_t minXForChunk = ii_MinX;
   
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      memset(ip_XBases[xIdx], 0x00, ii_XPerChunk * sizeof(uint32_t));
      
      // Set the first base for this chunk so we can use as an offset for the array in the GPU.
      ip_XBases[xIdx][0] = minXForChunk;
      minXForChunk += ii_XPerChunk;
   }
      
   // Reset these to 0
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
      ip_TermsInGroup[xIdx] = 0;
      
   term_t *terms = ip_LifchitzApp->GetTerms();

   idx = 0;

   while (terms[idx].x > 0)
   {
      xIdx = (terms[idx].x - ii_MinX) / ii_XPerChunk;
 
      minXForChunk = ip_XBases[xIdx][0];
 
      xOffset = terms[idx].x - minXForChunk;

      ip_XBases[xIdx][xOffset] = terms[idx].x;

      tIdx = ip_TermsInGroup[xIdx];

      ip_GpuTerms[xIdx][tIdx].x = terms[idx].x;
      ip_GpuTerms[xIdx][tIdx].y = terms[idx].y;
      ip_GpuSigns[xIdx][tIdx] = terms[idx].signs;

      ip_TermIndexes[xIdx][tIdx] = idx;
            
      ip_TermsInGroup[xIdx]++;
      
      idx++;
   }
   
   if (ii_XChunks > 1)
   {
      if (il_NextTermsBuild == 0)
         ip_LifchitzApp->WriteToConsole(COT_OTHER, "Building term groups");
      else
         ip_LifchitzApp->WriteToConsole(COT_OTHER, "Rebuilding term groups");
   }
   
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      // The last term for the group will have x = 0.
      minXForChunk = ip_XBases[xIdx][0];
      uint32_t maxXForChunk = minXForChunk + ii_XPerChunk;
         
      if (ii_XChunks > 1)
      {
         ip_LifchitzApp->WriteToConsole(COT_OTHER, "   Group %2u (%6u <= x < %6u) has %8u terms", 
            xIdx+1, minXForChunk, maxXForChunk, ip_TermsInGroup[xIdx]);
      }

      tIdx = ip_TermsInGroup[xIdx];
      
      // Make sure we have an "end of list"
      ip_GpuTerms[xIdx][tIdx].x = 0;
      ip_GpuTerms[xIdx][tIdx].y = 0;
      ip_GpuSigns[xIdx][tIdx] = 0;
   }
   
   if (ip_LifchitzApp->GetGpuWorkerCount() + ip_LifchitzApp->GetCpuWorkerCount() > 1)
      xfree(terms);;
}