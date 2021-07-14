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
#include "magiccl.h"
#include "pix_kernel.h"
#include <time.h>

#define   VECTOR_SIZE   2

PrimesInXGpuWorker::PrimesInXGpuWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   const char *magicSource[3];
   const char *source[10];
   
   ib_GpuWorker = true;
   
   ip_PrimesInXApp = (PrimesInXApp *) theApp;
   
   ii_MaxSteps = ip_PrimesInXApp->GetSteps();

   BuildMixedTerms();
   BuildSingleTerms();

   source[0] = pix_kernel;
   source[1] = 0;

   magicSource[0] = magic;
   magicSource[1] = 0;

   ip_MagicKernel = new Kernel(ip_PrimesInXApp->GetDevice(), "magic_kernel", magicSource);
   ip_PrimesInXKernel = new Kernel(ip_PrimesInXApp->GetDevice(), "pix_kernel", source);

   AllocatePrimeList(ip_PrimesInXKernel->GetWorkGroupSize());
   
   il_RemainderList = (int64_t *) xmalloc(ii_WorkSize*sizeof(int64_t));
   il_MagicNumber = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   il_MagicShift = (uint64_t *) xmalloc(ii_WorkSize*sizeof(uint64_t));
   ii_KernelDigitList = (uint32_t *) xmalloc(ii_MaxSteps*sizeof(uint32_t));
   il_FactorList = (uint64_t *) xmalloc(ii_MaxSteps*sizeof(uint64_t)); 

   ip_KAPrime     = new KernelArgument(ip_PrimesInXApp->GetDevice(), "prime", KA_HOST_TO_GPU, il_PrimeList, ii_WorkSize);
   ip_MKAMagicNumber = new KernelArgument(ip_PrimesInXApp->GetDevice(), "magic_number", KA_GPU_TO_HOST, il_MagicNumber, ii_WorkSize);
   ip_MKAMagicShift = new KernelArgument(ip_PrimesInXApp->GetDevice(), "magic_shift", KA_GPU_TO_HOST, il_MagicShift, ii_WorkSize);
   
   ip_MagicKernel->AddArgument(ip_KAPrime);
   ip_MagicKernel->AddArgument(ip_MKAMagicNumber);
   ip_MagicKernel->AddArgument(ip_MKAMagicShift);
   
   ip_KARemainders  = new KernelArgument(ip_PrimesInXApp->GetDevice(), "remainder", KA_BIDIRECTIONAL, il_RemainderList, ii_WorkSize);
   ip_SKAMagicNumber = new KernelArgument(ip_PrimesInXApp->GetDevice(), "magic_number", KA_HOST_TO_GPU, il_MagicNumber, ii_WorkSize);
   ip_SKAMagicShift = new KernelArgument(ip_PrimesInXApp->GetDevice(), "magic_shift", KA_HOST_TO_GPU, il_MagicShift, ii_WorkSize);
   ip_KADigits = new KernelArgument(ip_PrimesInXApp->GetDevice(), "digits", KA_HOST_TO_GPU, ii_KernelDigitList, ii_MaxSteps);  
   ip_KAFactors = new KernelArgument(ip_PrimesInXApp->GetDevice(), "factors", KA_BIDIRECTIONAL, il_FactorList, ii_MaxSteps);  

   ip_PrimesInXKernel->AddArgument(ip_KAPrime);
   ip_PrimesInXKernel->AddArgument(ip_SKAMagicNumber);
   ip_PrimesInXKernel->AddArgument(ip_SKAMagicShift);
   ip_PrimesInXKernel->AddArgument(ip_KARemainders);
   ip_PrimesInXKernel->AddArgument(ip_KADigits);
   ip_PrimesInXKernel->AddArgument(ip_KAFactors);
   
   ib_Initialized = true;
}

void  PrimesInXGpuWorker::CleanUp(void)
{
   delete ip_KAPrime;
   delete ip_MKAMagicNumber;
   delete ip_MKAMagicShift;
   
   delete ip_SKAMagicNumber;
   delete ip_SKAMagicShift;
   delete ip_KARemainders;
   delete ip_KAFactors;
   delete ip_KADigits;
   
   xfree(ii_KernelDigitList);
   xfree(il_RemainderList);
   xfree(il_FactorList);
   xfree(il_MagicNumber);
   xfree(il_MagicShift);
   
   for (uint32_t ii=0; ii<ii_MixedDigitListCount; ii++)
      xfree(ii_MixedDigitList[ii]);

   xfree(ii_MixedDigitList);

   for (uint32_t ii=0; ii<ii_SingleDigitListCount; ii++)
      xfree(ii_SingleDigitList[ii]);

   xfree(ii_SingleDigitList);
}

void  PrimesInXGpuWorker::TestMegaPrimeChunk(void)
{
   uint32_t ii, kk;
   int64_t  prime, rem;
   uint32_t term, plength;
   uint32_t digitListUsed;
   time_t   reportTime;
   uint32_t **digitList;
   
   // Once we reach 1G, then we can use larger terms to get
   // to the min term more quickly.
   if (il_PrimeList[0] > 1000000000)
   {
      digitList = ii_MixedDigitList;
      digitListUsed = ii_MixedDigitListUsed;
   }
   else
   {
      digitList = ii_SingleDigitList;
      digitListUsed = ii_SingleDigitListUsed;
   }
   
   ip_MagicKernel->Execute(ii_WorkSize);
   
   // Initialize the list of remainders
   for (ii=0; ii<ii_WorkSize; ii++)
      il_RemainderList[ii] = 0;
   
   reportTime = time(NULL) + 60;

   kk = 0;
   do
   {     
      for (ii=0; ii<ii_MaxSteps; ii++)
      {
         ii_KernelDigitList[ii] = digitList[kk][ii];
         il_FactorList[ii] = PMAX_MAX_62BIT;
      }

      ip_PrimesInXKernel->Execute(ii_WorkSize);

      // The first two entries will not have factors.
      for (ii=2; ii<ii_MaxSteps; ii++)
      {
         if (ii_KernelDigitList[ii] == ii_KernelDigitList[1])
            break;

         if (il_FactorList[ii] == PMAX_MAX_62BIT)
            continue;
         
         term = (ii - 1) + ii_KernelDigitList[0];

         prime = il_FactorList[ii];
         plength = 1;
         rem = prime;
         while (rem > 10)
         {
            plength++;
            rem /= 10;
         }
         
         if (plength == term)
            ip_PrimesInXApp->ReportPrime(prime, term);
         else
            if (ip_PrimesInXApp->ReportFactor(prime, term))
               VerifyFactor(prime, term);
      }

      kk++;

      if (kk < digitListUsed && ip_PrimesInXApp->IsInterrupted() && time(NULL) > reportTime)
      {
         ip_PrimesInXApp->WriteToConsole(COT_SIEVE, "Thread %d has completed %d of %d iterations", ii_MyId, kk, digitListUsed);
         reportTime = time(NULL) + 60;
      }
   } while (kk < digitListUsed);
   
   SetLargestPrimeTested(il_PrimeList[ii_WorkSize-1], ii_WorkSize);
}

void  PrimesInXGpuWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("PrimesInXGpuWorker::TestMiniPrimeChunk not implemented");
}

void  PrimesInXGpuWorker::BuildMixedTerms(void)
{
   uint32_t    ii, index;
   uint32_t    digitLength = 0;
   uint32_t    e9Terms, e1Terms;
   uint32_t    e9Term, e1Term;
   uint32_t   *digitList;

   e9Terms = (ip_PrimesInXApp->GetMinLength() / 9);
   e1Terms = ip_PrimesInXApp->GetMaxLength() - (9 * e9Terms);

   ii_MixedDigitListCount = 2 + (e9Terms + e1Terms) / (ii_MaxSteps - 4);

   ii_MixedDigitList = (uint32_t **) xmalloc(ii_MixedDigitListCount*sizeof(uint32_t *));
   
   for (ii=0; ii<ii_MixedDigitListCount; ii++)
      ii_MixedDigitList[ii] = (uint32_t *) xmalloc(2 + ii_MaxSteps*sizeof(uint32_t));

   index = 0;
   e1Term = e9Term = 0;

   // ii_MixedDigitList[index][0] is the total length of terms prior to this list
   ii_MixedDigitList[index][0] = 0;
   if (e9Terms == 0)
   {
      ii_MixedDigitList[index][1] = 10;
      e1Term = 2;
   }
   else
   {
      ii_MixedDigitList[index][1] = 1000000000;
      e9Term = 2;
   }

   digitList = ip_PrimesInXApp->Get1DigitTerms();
   for (ii=0; ii<ip_PrimesInXApp->GetMaxLength(); ii++)
   {
      if (ii < e9Terms * 9)
      {
         if (ii > 0 && ii % 9 == 0)
         {
            e9Term++;

            // If we don't have room for more terms, set the last term
            // to 1e9 to let the kernel know when to exit the loop, then
            // create the next digitList
            if (e9Term == (ii_MaxSteps - 2))
            {
               ii_MixedDigitList[index][e9Term] = 1000000000;
               index++;
               ii_MixedDigitList[index][0] = digitLength;
               ii_MixedDigitList[index][1] = 1000000000;
               e9Term = 2;
            }
         }

         ii_MixedDigitList[index][e9Term] *= 10;
         ii_MixedDigitList[index][e9Term] += digitList[ii];
         digitLength++;
      }
      else
      {
         if (ii > 0 && ii == e9Terms * 9)
         {
            e9Term++;
            ii_MixedDigitList[index][e9Term] = 1000000000;
            index++;
            ii_MixedDigitList[index][0] = digitLength;
            ii_MixedDigitList[index][1] = 10;
            e1Term = 2;
         }

         // If we don't have room for more terms, set the last term
         // to 1e1 to let the kernel know when to exit the loop, then
         // create the next digitList
         if (e1Term == (ii_MaxSteps - 2))
         {
            ii_MixedDigitList[index][e1Term] = 10;
            index++;
            ii_MixedDigitList[index][0] = digitLength;
            ii_MixedDigitList[index][1] = 10;
            e1Term = 2;
         }

         ii_MixedDigitList[index][e1Term] = digitList[ii];
         digitLength++;
         e1Term++;

      }
   }

   ii_MixedDigitList[index][e1Term] = 10;
   
   ii_MixedDigitListUsed = index + 1;
}

void  PrimesInXGpuWorker::BuildSingleTerms(void)
{
   uint32_t    ii, index;
   uint32_t    digitLength = 0;
   uint32_t    e1Terms;
   uint32_t    e1Term;
   uint32_t   *digitList;

   e1Terms = ip_PrimesInXApp->GetMaxLength();

   ii_SingleDigitListCount = 2 + e1Terms / (ii_MaxSteps - 4);

   ii_SingleDigitList = (uint32_t **) xmalloc(ii_SingleDigitListCount*sizeof(uint32_t *));
   
   for (ii=0; ii<ii_SingleDigitListCount; ii++)
      ii_SingleDigitList[ii] = (uint32_t *) xmalloc(2 + ii_MaxSteps*sizeof(uint32_t));

   index = 0;

   // ii_SingleDigitList[index][0] is the total length of terms prior to this list
   ii_SingleDigitList[index][0] = 0;
   ii_SingleDigitList[index][1] = 10;
   e1Term = 2;

   digitList = ip_PrimesInXApp->Get1DigitTerms();
   for (ii=0; ii<ip_PrimesInXApp->GetMaxLength(); ii++)
   {
      // If we don't have room for more terms, set the last term
      // to 1e1 to let the kernel know when to exit the loop, then
      // create the next digitList
      if (e1Term == (ii_MaxSteps - 2))
      {
         ii_SingleDigitList[index][e1Term] = 10;
         index++;
         ii_SingleDigitList[index][0] = digitLength;
         ii_SingleDigitList[index][1] = 10;
         e1Term = 2;
      }

      ii_SingleDigitList[index][e1Term] = digitList[ii];
      digitLength++;
      e1Term++;
   }

   ii_SingleDigitList[index][e1Term] = 10;
   
   ii_SingleDigitListUsed = index + 1;
}

void  PrimesInXGpuWorker::VerifyFactor(uint64_t prime, uint32_t n)
{
   // Note that p is limited to 2^52, so we are not using extended precision
   // in this function.

   // We know that a factor was found, but the call to pixsieve
   // doesn't tell us the term or the p for it, so we will use
   // the long way to find it.
   double   inverse, qd;
   uint32_t i;
   int64_t  p;
   int64_t  q, rem;
   
   uint32_t *terms = ip_PrimesInXApp->Get1DigitTerms();
      
   inverse = 1.0/prime;
   rem = 0;

   p = prime;
   for (i=0; i<n; i++)
   {
      qd = (((double)rem * 10.0) + (double)terms[i]);

      q = (int64_t) (qd * inverse);

      rem = (rem*10)+terms[i] - p*q;

      if (rem < 0)
         rem += p;
      else if (rem >= p)
         rem -= p;
   }

   if (rem != 0)
      FatalError("%" PRIu64" does not divide pix(%u) (remainder is %" PRIu64")", prime, n, rem);
}