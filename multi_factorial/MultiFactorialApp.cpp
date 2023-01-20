/* MultiFactorialApp.cpp -- (C) Mark Rodenkirch, October 2012

   This program finds factors of numbers in the form n!m+1 and n!m-1

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include "MultiFactorialApp.h"
#include "MultiFactorialWorker.h"
#include "../core/Parser.h"
#include "../core/MpArith.h"
#include "../sieve/primesieve.hpp"

#if defined(USE_OPENCL)
#include "MultiFactorialGpuWorker.h"
#define APP_NAME        "mfsievecl"
#elif defined(USE_METAL)
#include "MultiFactorialGpuWorker.h"
#define APP_NAME        "mfsievemtl"
#else
#define APP_NAME        "mfsieve"
#endif

#define APP_VERSION     "2.2"

#define BIT(n)          ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the CPUSieve framework for other applications.
App *get_app(void)
{
   return new MultiFactorialApp();
}

MultiFactorialApp::MultiFactorialApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of multi-factorials");
   SetLogFileName("mfsieve.log");

   ii_MultiFactorial = 1;
   ii_MinN = 0;
   ii_MaxN = 0;
   ii_CpuWorkSize = 50000;
   
   // We'll remove all even terms manually
   SetAppMinPrime(3);

   SetAppMaxPrime(PMAX_MAX_62BIT);
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuSteps = 100000;
   ii_MaxGpuFactors = 10;
#endif
}

MultiFactorialApp::~MultiFactorialApp()
{
   if (ip_Terms == NULL)
      return;

   for (uint32_t mf=0; mf<ii_MultiFactorial; mf++)
   {
      xfree(ip_Terms[mf].base);
      xfree(ip_Terms[mf].power);
   }
   
   xfree(ip_Terms);
}

void MultiFactorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=M           maximum n to search\n");
   printf("-m --multifactorial=m multifactorial, e.g. x!m where m = 3 --> x!!! (default 1)\n");

#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-S --step=S              max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactordensity=M  max number of factors to support per 1000 n per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  MultiFactorialApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:m:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');
   AppendLongOpt(longOpts, "multifactorial", required_argument, 0, 'm');

#if defined(USE_OPENCL) || defined(USE_METAL)
   shortOpts += "S:M:";
   
   AppendLongOpt(longOpts, "maxsteps",       required_argument, 0, 'S');
   AppendLongOpt(longOpts, "maxfactors",     required_argument, 0, 'M');
#endif
}

parse_t MultiFactorialApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'm':
         status = Parser::Parse(arg, 1, 1000000000, ii_MultiFactorial);
         break;
         
      case 'n':
         status = Parser::Parse(arg, 1, 1000000000, ii_MinN);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxN);
         break;

#if defined(USE_OPENCL) || defined(USE_METAL)
      case 'S':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxGpuSteps);
         break;
         
      case 'M':
         status = Parser::Parse(arg, 1, 1000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void MultiFactorialApp::ValidateOptions(void)
{
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
   
      iv_MinusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
      
      iv_PlusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);

      il_TermCount = 0;
      ProcessInputTermsFile(true);
   }
   else
   {        
      if (ii_MinN <= ii_MultiFactorial)
         FatalError("The value for -n must be greater than the value for -m");

      if (ii_MaxN <= ii_MinN)
         FatalError("The value for -N must be greater than the value for -n");

      iv_MinusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), true);
      
      iv_PlusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), true);
      
      il_TermCount = 2 * (ii_MaxN - ii_MinN + 1);

      // For even multi-factorials, remove all even terms
      if (ii_MultiFactorial % 2 == 0)
      {
         uint32_t n = ii_MinN;
         if (n % 2 == 0)
            n++;
         
         while (n <= ii_MaxN)
         {
            iv_PlusTerms[BIT(n)] = false;
            iv_MinusTerms[BIT(n)] = false;
            il_TermCount -= 2;
            
            n += 2;
         }
      }
   }

   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[50];
      
      if (ii_MultiFactorial == 1)
         snprintf(fileName, sizeof(fileName), "factorial.pfgw");
      else
         snprintf(fileName, sizeof(fileName), "mf_%d.pfgw", ii_MultiFactorial);
      
      is_OutputTermsFileName = fileName;
   }

   if (GetMinPrime() < GetMinN() && ii_MultiFactorial == 1)
   {
      WriteToConsole(COT_OTHER, "Setting min prime to %u since no primes below that will yield a factor", GetMinN());
      SetMinPrime(GetMinN());
   }
 
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuFactors = ii_MaxGpuFactors * GetGpuWorkGroups() * (ii_MaxN - ii_MinN) / 1000;
#endif

   FactorApp::ParentValidateOptions();

   BuildTerms();
   
   // The testing routine is optimized to test 4 primes at a time.
   while (ii_CpuWorkSize % 4 > 0)
      ii_CpuWorkSize++;
}

Worker *MultiFactorialApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#if defined(USE_OPENCL) || defined(USE_METAL)
   if (gpuWorker)
      return new MultiFactorialGpuWorker(id, this);
#endif
   
   return new MultiFactorialWorker(id, this);
}

void   MultiFactorialApp::BuildTerms(void)
{
   uint64_t *primes, thePrime;
   uint32_t  primeCount, mfIdx;
   primesieve::iterator   primeIterator;
   
   primeIterator.skipto(1, ii_MinN);
   
   primeCount = primesieve::count_primes(1, ii_MinN);
   primes = (uint64_t *) xmalloc(primeCount * sizeof(uint64_t));
   primeCount = 0;
   
   do
   {
      thePrime = primeIterator.next_prime();
      
      primes[primeCount] = thePrime;
      primeCount++;
   } while (thePrime < ii_MinN);
   
   ip_Terms = (terms_t *) xmalloc(ii_MultiFactorial * sizeof(terms_t));
   
   for (uint32_t mf=0; mf<ii_MultiFactorial; mf++)
   {
      ip_Terms[mf].base = (uint64_t *) xmalloc((primeCount+1) * sizeof(uint64_t));
      ip_Terms[mf].power = (uint32_t *) xmalloc((primeCount+1) * sizeof(uint32_t));
      
      memcpy(ip_Terms[mf].base, primes, primeCount * sizeof(uint64_t));
   }

   uint64_t *allTerms = (uint64_t *) xmalloc(ii_MinN * sizeof(uint64_t));
   
   for (uint32_t f=1; f<ii_MinN; f++)
      allTerms[f] = f;

   uint64_t powers[32];
   
   for (uint32_t pIdx=0; pIdx<primeCount; pIdx++)
   {
      uint32_t idx;
      
      // Build a table of powers for each prime from p^1 to p^x until p^x > minn
      powers[0] = 1;
      powers[1] = primes[pIdx];
      
      for (idx=2; powers[idx-2] * powers[1] < ii_MinN; idx++)
         powers[idx] = powers[idx-1] * powers[1];

      idx--;
            
      while (idx > 0)
      {
         uint64_t pIdxFactor = powers[idx];
         uint32_t tIdx = pIdxFactor;

         // We can stop when tIdx >= minn
         while (tIdx < ii_MinN)
         {      
            // If p^tIdx divides the term, then update the power by tIdx and divide the term by p^tIdx.
            // In other words for the term 2^11, we will have a factor of 2^10, but not of 2^8 since
            // we will have divided it by 2^10 already.
            if (allTerms[tIdx] % pIdxFactor == 0)
            {
               mfIdx = (tIdx % ii_MultiFactorial);
               
               ip_Terms[mfIdx].power[pIdx] += idx;
               allTerms[tIdx] /= pIdxFactor;
            }
            
            // If pIdxFactor = 27, then we go from 27 to 54 to 81, etc.
            // These are the only terms that could have 27 as a factor.
            tIdx += pIdxFactor;
         }
         
         // Decrease the power
         idx--;
      }
   }

   // At this point base[0]^power[0] * base{1]^power[1] * ... * base{primeCount]^power[primeCount] can
   // be used to compute (minn-1)!mf for each mf.  Now we can combine bases with the same power as long
   // as the product of the two bases does not exceed PMAX_MAX_62BIT.  I think we can actually go up to
   // 2^63 or even 2^64, but we'll choose to be safe because it won't have much of an impact on speed.
   for (mfIdx=0; mfIdx<ii_MultiFactorial; mfIdx++)
   {
      //printf("mf=%u\n", mfIdx);
      
      for (uint32_t i=0; i<primeCount; i++)
      {
         if (ip_Terms[mfIdx].power[i] == 0)
            continue;


      //   printf("%llu^%u * ", ip_Terms[mfIdx].base[i], ip_Terms[mfIdx].power[i]);
      }
      
      //printf("\n");
         
      for (uint32_t i=0; i<primeCount; i++)
      {
         if (ip_Terms[mfIdx].power[i] == 0)
            continue;
         
         for (uint32_t j=i+1; j<primeCount; j++)
         {
            if (ip_Terms[mfIdx].power[j])
               continue;
            
            if (ip_Terms[mfIdx].power[i] == ip_Terms[mfIdx].power[j])
            {               
               uint64_t mult = ip_Terms[mfIdx].base[i] * ip_Terms[mfIdx].base[j];
               
               // If we overflows then mult will be less than terms[mf].base[i]
               if ((double) ip_Terms[mfIdx].base[i] * (double) ip_Terms[mfIdx].base[j] < (double) PMAX_MAX_62BIT)
               {
                  ip_Terms[mfIdx].base[i] = mult;
                  ip_Terms[mfIdx].power[j] = 0;
               }
            }
         }         
      }

      // Now collapse the terms so that the last one has 
      for (uint32_t i=0; i<primeCount; i++)
      {
         if (ip_Terms[mfIdx].power[i] != 0)
            continue;

         bool filled = false;
         
         for (uint32_t j=i+1; j<primeCount; j++)
         {
            if (ip_Terms[mfIdx].power[j] > 0)
            {
               ip_Terms[mfIdx].base[i] = ip_Terms[mfIdx].base[j];
               ip_Terms[mfIdx].power[i] = ip_Terms[mfIdx].power[j];
               ip_Terms[mfIdx].power[j] = 0;
               filled = true;
               break;
            }
         }
         
         if (!filled)
            break;
      }
   }

   for (uint32_t mfIdx=0; mfIdx<ii_MultiFactorial; mfIdx++)
   {
      uint32_t count = 0;
      //printf("mf=%u\n", mfIdx);
      
      for (count=0; count<primeCount; count++)
      {
         if (ip_Terms[mfIdx].power[count] == 0)
            break;
         
         //printf("%llu^%u * ", ip_Terms[mfIdx].base[count], ip_Terms[mfIdx].power[count]);
      }
      
      ip_Terms[mfIdx].count = count;
      
      //printf("\n");
   }

   xfree(primes);
   xfree(allTerms);
}

void MultiFactorialApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000], *pos;
   uint32_t n;
   int32_t  c;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
  if (memcmp(buffer, "ABC $a#!$b", 9) && memcmp(buffer, "ABC $a!+$b", 10) &&
      sscanf(buffer, "ABC $a!%u$b", &ii_MultiFactorial) != 1)
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   pos = strstr(buffer, "//");
   if (pos)
      if (sscanf(pos+13, "%" SCNu64"", &sieveLimit) == 1)
         SetMinPrime(sieveLimit);

   if (!haveBitMap)
      ii_MinN = ii_MaxN = 0;
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
   
      if (sscanf(buffer, "%d %d", &n, &c) != 2)
         FatalError("Line %s is malformed", buffer);

      if (!ii_MaxN)
         ii_MinN = ii_MaxN = n;
            
      if (haveBitMap)
      {
         if (c == -1)
         {
            iv_MinusTerms[n - ii_MinN] = true;
            il_TermCount++;
         }
         
         if (c == +1)
         {
            iv_PlusTerms[n - ii_MinN] = true;
            il_TermCount++;
         }
      }
      else
      {
         if (ii_MinN > n) ii_MinN = n;
         if (ii_MaxN < n) ii_MaxN = n;
      }
   }

   fclose(fPtr);
}

bool MultiFactorialApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint32_t n, mf;
   int32_t  c;
   
   if (sscanf(term, "%u!%u%d", &n, &mf, &c) != 3)
      FatalError("Could not parse term %s", term);

   if (mf != ii_MultiFactorial)
      FatalError("Expected multifactorial of %u in factor, but found %u", ii_MultiFactorial, mf);
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
     
   VerifyFactor(theFactor, n, c);
   
   uint64_t bit = BIT(n);
   
   // No locking is needed because the Workers aren't running yet
   if (c == -1 && iv_MinusTerms[bit])
   {	
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      return true;
   }

   if (c == +1 && iv_PlusTerms[bit])
   {	
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      return true;
   }
      
   return false;
}

void MultiFactorialApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint64_t termsCounted = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC $a!%d$b // Sieved to %" PRIu64"\n", ii_MultiFactorial, largestPrime);

   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      if (iv_MinusTerms[n - ii_MinN])
      {
         fprintf(termsFile, "%d -1\n", n);
         termsCounted++;
      }

      if (iv_PlusTerms[n - ii_MinN])
      {
         fprintf(termsFile, "%d +1\n", n);
         termsCounted++;
      }
   }

   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void MultiFactorialApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   if (ii_MultiFactorial > 1)
      snprintf(extraText, maxTextLength, "%d <= n <= %d, multifactorial %u", ii_MinN, ii_MaxN, ii_MultiFactorial);
   else
      snprintf(extraText, maxTextLength, "%d <= n <= %d, factorial", ii_MinN, ii_MaxN);
}

bool MultiFactorialApp::ReportFactor(uint64_t theFactor, uint32_t n, int32_t c)
{
   uint32_t bit;
   bool     newFactor = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   VerifyFactor(theFactor, n, c);
   
   ip_FactorAppLock->Lock();

   bit = BIT(n);
   
   if (c == -1 && iv_MinusTerms[bit])
   {
      newFactor = true;
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(theFactor, "%u!%u-1", n, ii_MultiFactorial);
   }

   if (c == +1 && iv_PlusTerms[bit])
   {
      newFactor = true;
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(theFactor, "%u!%u+1", n, ii_MultiFactorial);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

void  MultiFactorialApp::VerifyFactor(uint64_t theFactor, uint32_t n, int32_t c)
{
   uint32_t currN = n;
   bool     termIsPrime = true;
   bool     isValid = false;
   char     buffer[100];
   
   uint64_t rem = 1, nextRem;
   
   while (termIsPrime && currN > 1)
   {
      nextRem = (rem * currN);
      
      // If rem * currN > UINT64_MAX, then it will truncate the result and
      // if truncated, then the term cannot be prime.
      if (nextRem < rem)
         termIsPrime = false;
      
      if (nextRem > theFactor + 1)
         termIsPrime = false;
      
      rem = nextRem;
      currN -= ii_MultiFactorial;
   }

   if (termIsPrime)
   {
      if (ii_MultiFactorial == 1)
         snprintf(buffer, sizeof(buffer), "%u!%+d is prime! (%" PRId64")", n, c, theFactor);
      else
         snprintf(buffer, sizeof(buffer), "%u!%u%+d is prime! (%" PRId64")", n, ii_MultiFactorial, c, theFactor);
         
      WriteToConsole(COT_OTHER, "%s", buffer);

      WriteToLog("%s", buffer);
      
      return;
   }

   if (n % ii_MultiFactorial == 0)
      currN = ii_MultiFactorial;
   else
      currN = (n % ii_MultiFactorial);
  
   MpArith  mp(theFactor);
   MpRes    pOne = mp.one();
   MpRes    mOne = mp.sub(mp.zero(), pOne);

   MpRes    resMf = mp.nToRes(ii_MultiFactorial);
   MpRes    ri = mp.nToRes(currN);
   MpRes    rf = ri;

   while (currN < n)
   {
      currN += ii_MultiFactorial;
      
      ri = mp.add(ri, resMf);
      rf = mp.mul(rf, ri);
   }

   if (c == -1)
      isValid = (rf == pOne);
      
   if (c == +1)
      isValid = (rf == mOne);
      
   if (isValid)
      return;
   
   if (ii_MultiFactorial == 1)
      snprintf(buffer, sizeof(buffer), "%" PRIu64" is not a factor of not a factor of %u!%+d", theFactor, n, c);
   else
      snprintf(buffer, sizeof(buffer), "%" PRIu64" is not a factor of not a factor of %u!%u%+d", theFactor, n, ii_MultiFactorial, c);

   FatalError("%s", buffer);
}
