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
   ii_MinY = ip_LifchitzApp->GetMinY();
   ii_MaxY = ip_LifchitzApp->GetMaxY();

   ii_MaxGpuFactors = ip_LifchitzApp->GetMaxGpuFactors();

   AllocateMemory();

   snprintf(defines[defineCount++], 50, "#define D_MAX_X_BASES %u", ii_XBasesPerGroup);
   snprintf(defines[defineCount++], 50, "#define D_MAX_Y_BASES %u", ii_YBasesPerGroup);
   snprintf(defines[defineCount++], 50, "#define D_MAX_FACTORS %u", ii_MaxGpuFactors);
      
   for (idx=0; idx<defineCount; idx++)
      preKernelSources[idx] = defines[idx];
   
   preKernelSources[idx] = 0;
   
   ip_Kernel = (GpuKernel *) ip_App->GetGpuDevice()->CreateKernel("lifchitz_kernel", lifchitz_kernel, preKernelSources);

   ip_LifchitzApp->SetGpuWorkGroupSize(ip_Kernel->GetWorkGroupSize());
   
   ii_PrimesInList = ip_LifchitzApp->GetGpuPrimesPerWorker();
   
   il_PrimeList    = (uint64_t *)  ip_Kernel->AddCpuArgument("primes", sizeof(uint64_t),  ii_PrimesInList);
   ii_KernelXBases = (uint32_t *)  ip_Kernel->AddCpuArgument("xbases", sizeof(uint32_t),  ii_XBasesPerGroup);
   ii_KernelYBases = (uint32_t *)  ip_Kernel->AddCpuArgument("ybases", sizeof(uint32_t),  ii_YBasesPerGroup);
   ip_TermList     = (gputerm_t *) ip_Kernel->AddCpuArgument("terms",  sizeof(gputerm_t), ii_TermsPerGroup);
   ii_FactorCount  = (uint32_t *)  ip_Kernel->AddSharedArgument("factorCount", sizeof(uint32_t), 10);
   ip_FactorList   = (factor_t *)  ip_Kernel->AddGpuArgument("factorList", sizeof(factor_t), ii_MaxGpuFactors);

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
   uint32_t idx, xIdx, yIdx;

   ii_XChunks = ip_LifchitzApp->GetXChunks();
   ii_YChunks = ip_LifchitzApp->GetYChunks();

   ii_XBasesPerGroup = 1 + ((1 + ii_MaxX - ii_MinX) / ii_XChunks);
   ii_YBasesPerGroup = 1 + ((1 + ii_MaxY - ii_MinY) / ii_YChunks);

   ip_XBases = (uint32_t **) xmalloc(ii_XChunks * sizeof(uint32_t *));
   ip_YBases = (uint32_t **) xmalloc(ii_YChunks * sizeof(uint32_t *));

   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
      ip_XBases[xIdx] = (uint32_t *) xmalloc(ii_XBasesPerGroup * sizeof(uint32_t));
   
   for (yIdx=0; yIdx<ii_YChunks; yIdx++)
      ip_YBases[yIdx] = (uint32_t *) xmalloc(ii_YBasesPerGroup * sizeof(uint32_t));

   ip_TermsInGroup = (uint32_t **) xmalloc(ii_XChunks * sizeof(uint32_t *));
   
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
      ip_TermsInGroup[xIdx] = (uint32_t *) xmalloc(ii_YChunks * sizeof(uint32_t));
   
   term_t *terms = ip_LifchitzApp->GetTerms();
   idx = 0;

   while (terms[idx].x > 0)
   {
      xIdx = (terms[idx].x - ii_MinX) / ii_XBasesPerGroup;
      yIdx = (terms[idx].y - ii_MinY) / ii_YBasesPerGroup;
      
      ip_TermsInGroup[xIdx][yIdx]++;

      idx++;
   }
   
   ii_TermsPerGroup = 0; 
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      for (yIdx=0; yIdx<ii_YChunks; yIdx++)
      {
         if (ip_TermsInGroup[xIdx][yIdx] > ii_TermsPerGroup)
            ii_TermsPerGroup = ip_TermsInGroup[xIdx][yIdx];
      }
   }
   
   ip_GpuTerms = (gputerm_t ***) xmalloc(ii_XChunks * sizeof(gputerm_t **));
   ip_TermIndexes = (uint64_t ***) xmalloc(ii_XChunks * sizeof(uint64_t **));

   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      ip_GpuTerms[xIdx] = (gputerm_t **) xmalloc(ii_YChunks * sizeof(gputerm_t *));
      ip_TermIndexes[xIdx] = (uint64_t **) xmalloc(ii_YChunks * sizeof(uint64_t *));

      for (yIdx=0; yIdx<ii_YChunks; yIdx++)
      {
         ip_GpuTerms[xIdx][yIdx] = (gputerm_t *) xmalloc((1 + ii_TermsPerGroup) * sizeof(gputerm_t));
         ip_TermIndexes[xIdx][yIdx] = (uint64_t *) xmalloc((1 + ii_TermsPerGroup) * sizeof(uint64_t));
      }
   }
}

void  LifchitzGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t  group = 0;
   uint32_t  xIdx, yIdx;
   time_t    reportTime = time(NULL) + 60;

   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      CreateTermGroups();
      
      il_NextTermsBuild = (il_PrimeList[ii_PrimesInList-1] << 1);
   }

   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
      for (yIdx=0; yIdx<ii_YChunks; yIdx++)
      {
         group++;
         
         if (ip_LifchitzApp->IsInterrupted() && time(NULL) > reportTime)
         {
            ip_LifchitzApp->WriteToConsole(COT_SIEVE, "Thread %d was interrupted.  Completed %u of %u iterations", ii_MyId, group, ii_XChunks * ii_YChunks);
            reportTime = time(NULL) + 60;
         }
          
         memcpy(ii_KernelXBases, ip_XBases[xIdx], ii_XBasesPerGroup * sizeof(uint32_t));
         memcpy(ii_KernelYBases, ip_YBases[xIdx], ii_YBasesPerGroup * sizeof(uint32_t));
         memcpy(ip_TermList,  ip_GpuTerms[xIdx][yIdx], ii_TermsPerGroup * sizeof(gputerm_t));

         ii_FactorCount[0] = 0;
         
         ip_Kernel->Execute(ii_PrimesInList);

         for (uint32_t fIdx=0; fIdx<ii_FactorCount[0]; fIdx++)
         {
            uint64_t tIdx =  ip_FactorList[fIdx].termIdx;
            uint32_t x = (uint32_t) ip_TermList[tIdx].x;
            uint32_t y = (uint32_t) ip_TermList[tIdx].y;
            int32_t  sign = (ip_FactorList[fIdx].sign == 1 ? +1 : -1);
            uint64_t theFactor = ip_FactorList[fIdx].factor;
            
            ip_LifchitzApp->ReportFactor(theFactor, x, y, sign, ip_TermIndexes[xIdx][yIdx][tIdx]);
            
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
   uint32_t idx, xIdx, yIdx, tIdx;
   uint32_t xOffset, yOffset;

   // These arrays will have zeros where we don't need to compute x^x or y^y
   uint32_t minXForChunk = ii_MinX;
   
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      memset(ip_XBases[xIdx], 0x00, ii_XBasesPerGroup * sizeof(uint32_t));
      
      // Set the first base for this chunk so we can use as an offset for the array in the GPU.
      ip_XBases[xIdx][0] = minXForChunk;
      minXForChunk += ii_XBasesPerGroup;
   }
      
   uint32_t minYForChunk = ii_MinY;
   for (yIdx=0; yIdx<ii_YChunks; yIdx++)
   {
      memset(ip_YBases[yIdx], 0x00, ii_YBasesPerGroup * sizeof(uint32_t));
         
      // Set the first base for minYForChunk chunk so we can use as an offset for the array in the GPU.
      ip_YBases[yIdx][0] = minYForChunk;
      minYForChunk += ii_YBasesPerGroup;
   }
   
   for (xIdx=0; xIdx<ii_XChunks; xIdx++)
   {
      // The last term for the group will have x = 0.
      for (yIdx=0; yIdx<ii_YChunks; yIdx++)
      {
         ip_TermsInGroup[xIdx][yIdx] = 0;
         memset(ip_GpuTerms[xIdx][yIdx], 0x00, (1 + ii_TermsPerGroup) * sizeof(gputerm_t));
      }
   }
   
   term_t *terms = ip_LifchitzApp->GetTerms();
   idx = 0;

   while (terms[idx].x > 0)
   {
      xIdx = (terms[idx].x - ii_MinX) / ii_XBasesPerGroup;
      yIdx = (terms[idx].y - ii_MinY) / ii_YBasesPerGroup;
 
      minXForChunk = ip_XBases[xIdx][0];
      minYForChunk = ip_YBases[yIdx][0];
 
      xOffset = terms[idx].x - minXForChunk;
      yOffset = terms[idx].y - minYForChunk;

      ip_XBases[xIdx][xOffset] = terms[idx].x;
      ip_YBases[yIdx][yOffset] = terms[idx].y;

      tIdx = ip_TermsInGroup[xIdx][yIdx];

      ip_GpuTerms[xIdx][yIdx][tIdx].x = terms[idx].x;
      ip_GpuTerms[xIdx][yIdx][tIdx].y = terms[idx].y;
      ip_TermIndexes[xIdx][yIdx][tIdx] = idx;
      
      if (terms[idx].sign == +1) ip_GpuTerms[xIdx][yIdx][tIdx].sign = 1;
      if (terms[idx].sign == -1) ip_GpuTerms[xIdx][yIdx][tIdx].sign = 2;
      
      ip_TermsInGroup[xIdx][yIdx]++;
      
      idx++;
   }
}