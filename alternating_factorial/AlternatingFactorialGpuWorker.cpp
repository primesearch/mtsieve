/* AlternatingFactorial.cpp -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "AlternatingFactorialGpuWorker.h"
#include "magiccl.h"
#include "af_kernel.h"

AlternatingFactorialGpuWorker::AlternatingFactorialGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   const char *magicSource[3];
   const char *source[10];
   
   ip_AlternatingFactorialApp = (AlternatingFactorialApp *) theApp;
   
   ib_GpuWorker = true;
   
   ii_MaxFactors = ip_AlternatingFactorialApp->GetStepN() + 10;
   
   source[0] = af_kernel;
   source[1] = 0;

   magicSource[0] = magic;
   magicSource[1] = 0;

   ip_MagicKernel = new Kernel(ip_AlternatingFactorialApp->GetDevice(), "magic_kernel", magicSource);

   ip_AlternatingFactorialKernel = new Kernel(ip_AlternatingFactorialApp->GetDevice(), "af_kernel", source);

   AllocatePrimeList(ip_AlternatingFactorialKernel->GetWorkGroupSize());

   il_RemainderList = (int64_t *) xmalloc(ii_WorkSize*sizeof(int64_t));
   il_NM1TermList = (int64_t *) xmalloc(ii_WorkSize*sizeof(int64_t));
   il_MagicNumber = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MagicShift = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_FactorList = (uint64_t *) xmalloc(ii_MaxFactors*sizeof(uint64_t));

   ip_KAPrime     = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_MKAMagicNumber = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "magic_number", KA_GPU_TO_HOST, il_MagicNumber, ii_WorkSize);
   ip_MKAMagicShift = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "magic_shift", KA_GPU_TO_HOST, il_MagicShift, ii_WorkSize);

   ip_FKAMagicNumber = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "magic_number", KA_HOST_TO_GPU, il_MagicNumber, ii_WorkSize);
   ip_FKAMagicShift = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "magic_shift", KA_HOST_TO_GPU, il_MagicShift, ii_WorkSize);
   ip_KARemainder = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "remainder", KA_BIDIRECTIONAL, il_RemainderList, ii_WorkSize);
   ip_KANM1Term   = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "nm1_term", KA_BIDIRECTIONAL, il_NM1TermList, ii_WorkSize);
   ip_KANRange    = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "n_range", KA_HOST_TO_GPU, (int32_t *) &ii_NRange, 10);
   ip_KAFactor    = new KernelArgument(ip_AlternatingFactorialApp->GetDevice(), "factors", KA_BIDIRECTIONAL, il_FactorList, ii_MaxFactors);

   ip_MagicKernel->AddArgument(ip_KAPrime);
   ip_MagicKernel->AddArgument(ip_MKAMagicNumber);
   ip_MagicKernel->AddArgument(ip_MKAMagicShift);

   ip_AlternatingFactorialKernel->AddArgument(ip_KAPrime);
   ip_AlternatingFactorialKernel->AddArgument(ip_FKAMagicNumber);
   ip_AlternatingFactorialKernel->AddArgument(ip_FKAMagicShift);
   ip_AlternatingFactorialKernel->AddArgument(ip_KARemainder);
   ip_AlternatingFactorialKernel->AddArgument(ip_KANM1Term);
   ip_AlternatingFactorialKernel->AddArgument(ip_KANRange);
   ip_AlternatingFactorialKernel->AddArgument(ip_KAFactor);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  AlternatingFactorialGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_MKAMagicNumber;
   delete ip_MKAMagicShift;
   
   delete ip_FKAMagicNumber;
   delete ip_FKAMagicShift;
   delete ip_KARemainder;
   delete ip_KANM1Term;
   delete ip_KANRange;
   delete ip_KAFactor;

   xfree(il_RemainderList);
   xfree(il_NM1TermList);
   xfree(il_FactorList);
   xfree(il_MagicNumber);
   xfree(il_MagicShift);
}

void  AlternatingFactorialGpuWorker::TestMegaPrimeChunk(void)
{   
   uint32_t ii, kk, kkMax;
   int32_t  maxN = ip_AlternatingFactorialApp->GetMaxN();
   uint32_t stepN = ip_AlternatingFactorialApp->GetStepN();
   int64_t  prime;
   uint32_t term;
   time_t   reportTime;
   
   ip_MagicKernel->Execute(ii_WorkSize);
   ii_NRange[0] = 2;

   // Initialize the list of remainders
   for (ii=0; ii<ii_WorkSize; ii++)
   {
      il_RemainderList[ii] = 1;
      il_NM1TermList[ii] = 1;
   }

   // Note the use of < rather than <= both here and in the kernel
   maxN += 1;

   reportTime = time(NULL) + 60;

   kk = 0;
   kkMax = 1 + ((maxN - 2) / stepN);
   
   do
   {
      ii_NRange[1] = ii_NRange[0] + stepN;

      if (ii_NRange[1] > maxN)
         ii_NRange[1] = maxN;
      
      for (ii=0; ii<ii_MaxFactors; ii++)
         il_FactorList[ii] = PMAX_MAX_62BIT;

      ip_AlternatingFactorialKernel->Execute(ii_WorkSize);

      for (ii=0; ii<stepN; ii++)
      {
         if (il_FactorList[ii] == PMAX_MAX_62BIT)
            continue;
         
         term = ii + ii_NRange[0];
         prime = il_FactorList[ii];

         ip_AlternatingFactorialApp->ReportFactor(prime, term);
      }

      ii_NRange[0] = ii_NRange[1];

      ii_NRange[1] = ii_NRange[0] + stepN;

      if (ii_NRange[1] > maxN)
         ii_NRange[1] = maxN;

      ii_NRange[2] = 0;

      kk++;
      if (ii_NRange[0] < maxN && ip_AlternatingFactorialApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_AlternatingFactorialApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, kk, kkMax);
         reportTime = time(NULL) + 60;
      }
   } while (ii_NRange[0] < maxN);

   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  AlternatingFactorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("AlternatingFactorialGpuWorker::TestMiniPrimeChunk not implemented");
}

void AlternatingFactorialGpuWorker::ValidateFactor(int64_t p, uint32_t term)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.
   
   double   inverse, qd;
   uint32_t n;
   int64_t  up = p;
   int64_t  q, rem;
   int64_t  nm1term;
   
   inverse = 1.0 / p;
   nm1term = 1;
   rem = 1;
   
   for (n=2; n<=term; n++)
   {
      qd = ((double) rem * (double) n);

      q = (int64_t) (qd * inverse);

      rem = (rem*n) - p*q;

      if (rem < 0)
         rem += p;
      else if (rem >= up)
         rem -= p;
      
      // Set nm1term to af(n) % p
      if (rem > nm1term)
         nm1term = rem - nm1term;
      else
         nm1term = (rem + p - nm1term);
   }

   if (rem != 0)
      FatalError("%" PRIu64" does not divide af(n) (remainder is %" PRIu64")", p, term, rem);

}

