/* FactorialSieve.cpp -- (C) Mark Rodenkirch, September 2016

   This class sets up the call to the FSieve GPU function and parses the output
   from the GPU to determine if we have a factor.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "MultiFactorialGpuWorker.h"
#include "magiccl.h"
#include "mf_kernel.h"
#include <time.h>

MultiFactorialGpuWorker::MultiFactorialGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   char        define1[40];
   char        define2[40];
   const char *magicSource[3];
   const char *source[10];

   ip_MultiFactorialApp = (MultiFactorialApp *) theApp;

   ii_StepN = ip_MultiFactorialApp->GetStepN();
   ii_MultiFactorial = ip_MultiFactorialApp->GetMultiFactorial();
   
   ib_GpuWorker = true;

   sprintf(define1, "#define D_ADDER %d\n", ii_MultiFactorial);
   
   // Change if we add support to the CPU worker
   sprintf(define2, "\n");
   //sprintf(define2, "#define D_INVERSE\n");
      

   source[0] = define1;
   source[1] = define2;
   source[2] = mf_kernel;
   source[3] = 0;

   magicSource[0] = magic;
   magicSource[1] = 0;

   
   ip_MagicKernel = new Kernel(ip_MultiFactorialApp->GetDevice(), "magic_kernel", magicSource);

   ip_FactorialKernel = new Kernel(ip_MultiFactorialApp->GetDevice(), "mf_kernel", source);

   AllocatePrimeList(ip_FactorialKernel->GetWorkGroupSize());

   il_RemainderList = (int64_t *) xmalloc(ii_WorkSize*sizeof(int64_t));
   il_MagicNumber = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MagicShift = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MinusFactorList = (uint64_t *) xmalloc(ii_StepN*sizeof(uint64_t));
   il_PlusFactorList = (uint64_t *) xmalloc(ii_StepN*sizeof(uint64_t));

   ip_KAPrime        = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_MKAMagicNumber = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "magic_number", KA_GPU_TO_HOST, il_MagicNumber, ii_WorkSize);
   ip_MKAMagicShift  = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "magic_shift", KA_GPU_TO_HOST, il_MagicShift, ii_WorkSize);

   ip_FKAMagicNumber = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "magic_number", KA_HOST_TO_GPU, il_MagicNumber, ii_WorkSize);
   ip_FKAMagicShift  = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "magic_shit", KA_HOST_TO_GPU, il_MagicShift, ii_WorkSize);
   ip_KARemainder    = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "remainder", KA_BIDIRECTIONAL, il_RemainderList, ii_WorkSize);
   ip_KANRange       = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "n_range", KA_HOST_TO_GPU, (int32_t *) &ii_NRange, 10);
   ip_KAMinusFactor  = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "factor", KA_BIDIRECTIONAL, il_MinusFactorList, ii_StepN);
   ip_KAPlusFactor   = new KernelArgument(ip_MultiFactorialApp->GetDevice(), "factor", KA_BIDIRECTIONAL, il_PlusFactorList, ii_StepN);

   ip_MagicKernel->AddArgument(ip_KAPrime);
   ip_MagicKernel->AddArgument(ip_MKAMagicNumber);
   ip_MagicKernel->AddArgument(ip_MKAMagicShift);

   ip_FactorialKernel->AddArgument(ip_KAPrime);
   ip_FactorialKernel->AddArgument(ip_FKAMagicNumber);
   ip_FactorialKernel->AddArgument(ip_FKAMagicShift);
   ip_FactorialKernel->AddArgument(ip_KARemainder);
   ip_FactorialKernel->AddArgument(ip_KANRange);
   ip_FactorialKernel->AddArgument(ip_KAMinusFactor);
   ip_FactorialKernel->AddArgument(ip_KAPlusFactor);
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  MultiFactorialGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_KARemainder;
   delete ip_MKAMagicNumber;
   delete ip_MKAMagicShift;
   delete ip_FKAMagicNumber;
   delete ip_FKAMagicShift;
   delete ip_KANRange;
   delete ip_KAMinusFactor;
   delete ip_KAPlusFactor;

   xfree(il_MinusFactorList);
   xfree(il_PlusFactorList);
   xfree(il_MagicNumber);
   xfree(il_MagicShift);
   xfree(il_RemainderList);
}

void  MultiFactorialGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii;
   uint32_t maxN = ip_MultiFactorialApp->GetMaxN();
   int64_t  il;
   uint32_t n_seed, term;
   int32_t  iteration = 0, maxIterations;
   time_t   reportTime;

   ip_MagicKernel->Execute(ii_WorkSize);

   // This is an approximation
   maxIterations = ii_MultiFactorial * (1 + maxN) / ii_StepN;

   // For normal factorials, this loop will be iterated one.
   // For multi-factorials, it will be iterated for the value of the multi-factorial.  So
   // if this the multi-factorial is 3, i.e. n!3, then this is how it iterates:
   //    n_seed=1 will evaluate 1!3, 4!3, 7!3, etc.  (7!3 = 7*4*1)
   //    n_seed=2 will evaluate 2!3, 5!3, 8!3, etc.  (8!3 = 8*5*2)
   //    n_seed=3 will evaluate 3!3, 6!3, 9!3, etc.  (9!3 = 9*6*3)
   // For inverse multi-factorials, it will be similar.  If the inverse multi-factorial
   // is 3, i.e. n!/n!3, then this is how it iterates:
   //    n_seed=1 will evaluate 1!/1!3, 4!/4!3, 7!/7!3, etc.  (7!/7!3 = 6*5*3*2, i.e not multiplying by 7, 4, and 1)
   //    n_seed=2 will evaluate 2!/2!3, 5!/5!3, 8!/8!3, etc.  (8!/8!3 = 7*6*4*3*1, i.e. not multiplying by 8, 5, and 2)
   //    n_seed=3 will evaluate 3!/3!3, 6!/6!3, 9!/9!3, etc.  (9!/9!3 = 8*7*5*4*2*1, i.e. not multiplying by 9, 6, and 3))
   for (n_seed=1; n_seed<=ii_MultiFactorial; n_seed++)
   {
      ii_NRange[0] = n_seed;

      // Initialize the list of remainders
      for (il=0; il<ii_WorkSize; il++)
         il_RemainderList[il] = 1;

      // Note the use of < rather than <= both here and in the kernel
      maxN += 1;

      reportTime = time(NULL) + 60;

      while (ii_NRange[0] < maxN)
      {
         for (ii=0; ii<ii_StepN; ii++)
         {
            il_MinusFactorList[ii] = PMAX_MAX_62BIT;
            il_PlusFactorList[ii] = PMAX_MAX_62BIT;
         }

         ii_NRange[1] = ii_NRange[0] + (ii_StepN * ii_MultiFactorial);
         
         if (ii_NRange[1] > maxN)
            ii_NRange[1] = maxN;

         ip_FactorialKernel->Execute(ii_WorkSize);

         for (ii=0; ii<ii_StepN; ii++)
         {
            if (il_MinusFactorList[ii] != PMAX_MAX_62BIT)
            {
               term = ii*ii_MultiFactorial + ii_NRange[0];
               CheckForPrimeOrFactor(il_MinusFactorList[ii], term, -1);
            }
            
            if (il_PlusFactorList[ii] != PMAX_MAX_62BIT)
            {
               term = ii*ii_MultiFactorial + ii_NRange[0];
               CheckForPrimeOrFactor(il_PlusFactorList[ii], term, +1);
            }
         }

         ii_NRange[0] = ii_NRange[1];

         iteration++;

         if (iteration < maxIterations && ip_MultiFactorialApp->IsInterrupted() && time(NULL) > reportTime)
         {
            ip_MultiFactorialApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, iteration, maxIterations);
            reportTime = time(NULL) + 60;
         }  
      }
   }
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  MultiFactorialGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("MultiFactorialGpuWorker::TestMiniPrimeChunk not implemented");
}

void  MultiFactorialGpuWorker::CheckForPrimeOrFactor(int64_t p, uint32_t n, int32_t c)
{
   uint32_t term;
   int64_t  value;
   long double pd, vd;

   term = n;
   value = 1;
   
   pd = (long double) (p + 1);
   vd = (long double) value;   
   
   while (term > 0)
   {
      vd *= (long double) term;

      // If the current value is greater than the
      // the max, then the value cannot be prime.      
      if (vd > pd)
      {
         value = 0;
         break;
      }

      value *= term;
      term -= ii_MultiFactorial;
   }

   if (value == p+1)
   {
      ip_MultiFactorialApp->ReportPrime(p, n, -1);
      return;
   }
      
   if (value == p-1)
   {
      ip_MultiFactorialApp->ReportPrime(p, n, +1);
      return;
   }
                  
   if (ip_MultiFactorialApp->ReportFactor(p, n, c))
      ValidateFactor(p, n, c);
}

void  MultiFactorialGpuWorker::ValidateFactor(int64_t p, uint32_t n, int32_t c)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.

   double   inverse, qd;
   int32_t  term;
   int64_t  q, rem;

   inverse = 1.0/p;
   rem = 1;

   term = n;
   
   while (term > 1)
   {
      qd = ((double) rem * (double) term);

      q = (int64_t) (qd * inverse);

      rem = (rem*term) - p*q;

      if (rem < 0)
         rem += p;
      else if (rem >= (int64_t) p)
         rem -= p;
      
      term -= ii_MultiFactorial;
   }

   if (c == -1 && rem != +1)
   {
      if (ii_MultiFactorial == 1)
         FatalError("%" PRIu64" does not divide %u!-1 (remainder is %" PRIu64")", p, n, rem);
      else
         FatalError("%" PRIu64" does not divide %u!%u-1 (remainder is %" PRIu64")", p, n, ii_MultiFactorial, rem);
   }
   
   if (c == +1 && rem != p-1)
   {
      if (ii_MultiFactorial == 1)
         FatalError("%" PRIu64" does not divide %u!+1 (remainder is %" PRIu64")", p, n, rem);
      else
         FatalError("%" PRIu64" does not divide %u!%u+1 (remainder is %" PRIu64")", p, n, ii_MultiFactorial, rem);
   }
}

