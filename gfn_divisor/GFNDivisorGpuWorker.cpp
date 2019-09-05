/* GFNDivisorWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>

#include <cinttypes>
#include "GFNDivisorGpuWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "gfn_kernel.h"

GFNDivisorGpuWorker::GFNDivisorGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{    
   const char *source[10];
   char        define1[50];
   char        define2[50];
   char        define3[50];
   char        define4[50];
   char        define5[50];
   
   ip_GFNDivisorApp = (GFNDivisorApp *) theApp;
   
   ib_GpuWorker = true;
   
   il_KCount = ip_GFNDivisorApp->GetKCount();
   ii_NCount = ip_GFNDivisorApp->GetNCount();
   
   il_MinK = ip_GFNDivisorApp->GetMinK();
   il_MaxK = ip_GFNDivisorApp->GetMaxK();
   
   ii_MinN = ip_GFNDivisorApp->GetMinN();
   ii_MaxN = ip_GFNDivisorApp->GetMaxN();
   
   sprintf(define1, "#define N_MIN %u\n", ii_MinN);
   sprintf(define2, "#define N_MAX %u\n", ii_MaxN);
   sprintf(define3, "#define N_COUNT %u\n", ii_NCount);
   sprintf(define4, "#define K_MIN %" PRIu64"\n", il_MinK);
   sprintf(define5, "#define K_MAX %" PRIu64"\n", il_MaxK);
   
   source[0] = define1;
   source[1] = define2;
   source[2] = define3;
   source[3] = define4;
   source[4] = define5;
   source[5] = gfn_kernel;
   source[6] = 0;

   ip_GFNDivisorKernel = new Kernel(ip_GFNDivisorApp->GetDevice(), "gfn_kernel", source);

   AllocatePrimeList(ip_GFNDivisorKernel->GetWorkGroupSize());

   il_KList = (uint64_t *) xmalloc(ii_WorkSize*ii_NCount*sizeof(int64_t));

   ip_KAPrime     = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_KAKList     = new KernelArgument(ip_GFNDivisorApp->GetDevice(), "ks", KA_GPU_TO_HOST, il_KList, ii_WorkSize*ip_GFNDivisorApp->GetNCount());

   ip_GFNDivisorKernel->AddArgument(ip_KAPrime);
   ip_GFNDivisorKernel->AddArgument(ip_KAKList);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  GFNDivisorGpuWorker::CleanUp(void)
{  
   delete ip_KAPrime;
   delete ip_KAKList;

   xfree(il_KList);
}

void  GFNDivisorGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t index, pIndex, nIndex;
   
   ip_GFNDivisorKernel->Execute(ii_WorkSize);

   for (pIndex=0; pIndex<ii_WorkSize; pIndex++)
   {
      index = pIndex * ii_NCount;
      
      for (nIndex=0; nIndex<ii_NCount; nIndex++)
      {
         uint64_t k = il_KList[index + nIndex];
         
         if (k > 0)
         {
            uint32_t n = nIndex + ii_MinN;
            uint64_t prime = il_PrimeList[pIndex];
            
            if (ip_GFNDivisorApp->ReportFactor(prime, k, n))
               VerifyFactor(k, n, prime);
         }
      }
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  GFNDivisorGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("GFNDivisorGpuWorker::TestMiniPrimeChunk not implemented");
}

void  GFNDivisorGpuWorker::VerifyFactor(uint64_t k, uint32_t n, uint64_t prime)
{
   uint64_t rem, bPowN;

   fpu_push_1divp(prime);

   bPowN = fpu_powmod(2, n, prime);
   
   rem = fpu_mulmod(bPowN, k, prime);

   fpu_pop();
   
   if (rem != prime - 1)
      FatalError("%" PRIu64"*2^%u+1 mod %" PRIu64" = %" PRIu64"", k, n, prime, rem + 1);
}

