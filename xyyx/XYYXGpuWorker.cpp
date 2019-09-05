/* XYYXGpuWorker.cpp -- (C) Mark Rodenkirch, September 2012

   This class sets up the call to the XYYXGpuWorker GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>

#include "XYYXGpuWorker.h"
#include "magiccl.h"
#include "xyyx_kernel.h"
#include "../x86_asm/fpu-asm-x86.h"

XYYXGpuWorker::XYYXGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   const char *magicSource[3];
   const char *xyyxSource[5];
   
   ib_GpuWorker = true;
   
   ip_XYYXApp = (XYYXApp *) theApp;

   magicSource[0] = magic;
   magicSource[1] = 0;
   
   xyyxSource[0] = xyyx_kernel;
   xyyxSource[1] = 0;

   ip_MagicKernel = new Kernel(ip_XYYXApp->GetDevice(), "magic_kernel", magicSource);

   ip_XYYXKernel = new Kernel(ip_XYYXApp->GetDevice(), "xyyx_kernel", xyyxSource);

   AllocatePrimeList(ip_XYYXKernel->GetWorkGroupSize());
   
   // If there are more y for any x than the number of steps, then we'll give an error
   // to the user at runtime when building the term groups.
   ii_MinGroupSize = ip_XYYXApp->GetGpuSteps();
   ii_MaxGroupSize = 2 * ii_MinGroupSize;
   
   ii_MaxTermCount = ip_XYYXApp->GetTermCount() + 2 * ip_XYYXApp->GetXCount() + 2 * ip_XYYXApp->GetYCount();
   
   // Ensure the buffer big enough to avoid issues when executing the kernel
   if (ii_MaxTermCount < ii_MaxGroupSize)
      ii_MaxTermCount = ii_MaxGroupSize;

   // Allocate enough memory to hold all of the terms.
   ii_Terms = (uint32_t *) xmalloc(ii_MaxTermCount*sizeof(int32_t));
   
   il_PrimeList = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MagicNumber = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MagicShift = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MinusFactorList = (uint64_t *) xmalloc(ii_MaxGroupSize*sizeof(uint64_t));
   il_PlusFactorList = (uint64_t *) xmalloc(ii_MaxGroupSize*sizeof(uint64_t));
   
   ip_KAPrime        = new KernelArgument(ip_XYYXApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_MKAMagicNumber = new KernelArgument(ip_XYYXApp->GetDevice(), "magic_number", KA_GPU_TO_HOST, il_MagicNumber, ii_WorkSize);
   ip_MKAMagicShift  = new KernelArgument(ip_XYYXApp->GetDevice(), "magic_shift", KA_GPU_TO_HOST, il_MagicShift, ii_WorkSize);

   ip_XYYXKAMagicNumber = new KernelArgument(ip_XYYXApp->GetDevice(), "magic_number", KA_HOST_TO_GPU, il_MagicNumber, ii_WorkSize);
   ip_XYYXKAMagicShift  = new KernelArgument(ip_XYYXApp->GetDevice(), "magic_shift", KA_HOST_TO_GPU, il_MagicShift, ii_WorkSize);

   // We will change the pointer before executing the kernel
   ip_KATerms        = new KernelArgument(ip_XYYXApp->GetDevice(), "terms", KA_HOST_TO_GPU, ii_Terms, ii_MaxGroupSize);
   ip_KAMinusFactors = new KernelArgument(ip_XYYXApp->GetDevice(), "minus_factors", KA_BIDIRECTIONAL, il_MinusFactorList, ii_MaxGroupSize);
   ip_KAPlusFactors  = new KernelArgument(ip_XYYXApp->GetDevice(), "plus_factors", KA_BIDIRECTIONAL, il_PlusFactorList, ii_MaxGroupSize);

   ip_MagicKernel->AddArgument(ip_KAPrime);
   ip_MagicKernel->AddArgument(ip_MKAMagicNumber);
   ip_MagicKernel->AddArgument(ip_MKAMagicShift);

   ip_XYYXKernel->AddArgument(ip_KAPrime);
   ip_XYYXKernel->AddArgument(ip_XYYXKAMagicNumber);
   ip_XYYXKernel->AddArgument(ip_XYYXKAMagicShift);
   ip_XYYXKernel->AddArgument(ip_KATerms);
   ip_XYYXKernel->AddArgument(ip_KAMinusFactors);
   ip_XYYXKernel->AddArgument(ip_KAPlusFactors);

   // The thread can't start until initialization is done
   ib_Initialized = true;

   il_NextTermsBuild = 0;
}

void  XYYXGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;

   delete ip_MKAMagicNumber;
   delete ip_MKAMagicShift;
   
   delete ip_XYYXKAMagicNumber;
   delete ip_XYYXKAMagicShift;

   delete ip_KATerms;
   delete ip_KAPlusFactors;
   delete ip_KAMinusFactors;

   delete ip_MagicKernel;
   delete ip_XYYXKernel;

   xfree(ii_Terms);
   xfree(il_MagicNumber);
   xfree(il_MagicShift);
   xfree(il_MinusFactorList);
   xfree(il_PlusFactorList);
}

void  XYYXGpuWorker::TestMegaPrimeChunk(void)
{
   uint64_t  factor;
   uint32_t  group;
   uint32_t  ii, idx;
   uint32_t  x, y;
   time_t    reportTime = time(NULL) + 60;

   // Every once in a while rebuild the term lists as it will have fewer entries
   // which will speed up testing for the next range of p.
   if (il_PrimeList[0] > il_NextTermsBuild)
   {
      ii_Groups = ip_XYYXApp->GetTerms(ii_Terms, ii_MinGroupSize, ii_MaxGroupSize);
      
      il_NextTermsBuild = (il_PrimeList[0] << 1);
   }
      
   // Compute magic numbers
   ip_MagicKernel->Execute(ii_WorkSize);

   for (group=0; group<ii_Groups; group++) {
      if (ip_XYYXApp->IsInterrupted() || time(NULL) > reportTime) {
         ip_XYYXApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, group, ii_Groups);
         reportTime = time(NULL) + 60;
      }
      
      for (ii=0; ii<ii_MaxGroupSize; ii++)
      {
         il_MinusFactorList[ii] = PMAX_MAX_62BIT;
         il_PlusFactorList[ii] = PMAX_MAX_62BIT;
      }
      
      idx = group * ii_MaxGroupSize;
      ip_KATerms->SetHostMemory(&ii_Terms[idx]);

      ip_XYYXKernel->Execute(ii_WorkSize);

      ii = 0;
      while (ii_Terms[idx] != 0)
      {
         x = ii_Terms[idx];         
         idx++;
         ii++;
         
         while (ii_Terms[idx] != 0)
         {
            y = ii_Terms[idx];
            
            if (il_MinusFactorList[ii] != PMAX_MAX_62BIT)
            {
               factor = il_MinusFactorList[ii];
               
               if (ip_XYYXApp->ReportFactor(factor, x, y, -1))
                  VerifyFactor(factor, x, y, -1);
            }
                        
            if (il_PlusFactorList[ii] != PMAX_MAX_62BIT)
            {
               factor = il_PlusFactorList[ii];
               
               if (ip_XYYXApp->ReportFactor(factor, x, y, +1))
                  VerifyFactor(factor, x, y, +1);
            }
            
            // Go to the next term
            idx++;
            ii++;
         }
         
         // Skip the 0 (marking the end of y for this x
         idx++;
         ii++;
      }
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  XYYXGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("XYYXGpuWorker::TestMiniPrimeChunk not implemented");
}

void  XYYXGpuWorker::VerifyFactor(uint64_t p, uint32_t x, uint32_t y, int32_t c)
{
   uint64_t xPowY, yPowX;
   
   fpu_push_1divp(p);
      
   xPowY = fpu_powmod(x, y, p);
   yPowX = fpu_powmod(y, x, p);
   
   fpu_pop();
   
   if (c == -1 && xPowY != yPowX)
      FatalError("%" PRIu64" does not divide %u^%u-%u^%u", p, x, y, y, x);
   
   if (c == +1 && xPowY + yPowX != p)
      FatalError("%" PRIu64" does not divide %u^%u+%u^%u", p, x, y, y, x);
}