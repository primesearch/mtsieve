/* PrimorialApp.cpp -- (C) Mark Rodenkirch, July 2018

   PrimorialWorker/Wall-Sun-Sun Search OpenCL application

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include "PrimorialApp.h"
#include "PrimorialWorker.h"
#ifdef HAVE_GPU_WORKERS
#include "PrimorialGpuWorker.h"
#endif
#include "../core/Parser.h"
#include "../sieve/primesieve.hpp"
#include "../x86_asm/fpu-asm-x86.h"

#define APP_NAME        "psieve"
#define APP_VERSION     "1.4"

#define BIT(primorial)  ((primorial) - ii_MinPrimorial)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the CPUSieve framework for other applications.
App *get_app(void)
{
   return new PrimorialApp();
}

PrimorialApp::PrimorialApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of primorials");
   SetLogFileName("psieve.log");

   ii_MinPrimorial = 100;
   ii_MaxPrimorial = 0;
   ii_CpuWorkSize = 50000;
   
   // No reason to support smaller primorials since they are all known
   SetAppMinPrime(100);
   
   // This is because the assembly code is using SSE to do the mulmods
   SetAppMaxPrime(PMAX_MAX_52BIT);
   
#ifdef HAVE_GPU_WORKERS
   ii_MaxGpuSteps = 50000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}

void PrimorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum primorial to search\n");
   printf("-N --maxn=M           maximum primorial to search\n");

#ifdef HAVE_GPU_WORKERS
   printf("-S --step=S           max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  PrimorialApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');

#ifdef HAVE_GPU_WORKERS
   shortOpts += "S:M:";
   
   AppendLongOpt(longOpts, "maxsteps",       required_argument, 0, 'S');
   AppendLongOpt(longOpts, "maxfactors",     required_argument, 0, 'M');
#endif
}

parse_t PrimorialApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {         
      case 'n':
         status = Parser::Parse(arg, 100, 1000000000, ii_MinPrimorial);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 100, 1000000000, ii_MaxPrimorial);
         break;

#ifdef HAVE_GPU_WORKERS
      case 'S':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxGpuSteps);
         break;
         
      case 'M':
         status = Parser::Parse(arg, 10, 1000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void PrimorialApp::ValidateOptions(void)
{
   uint32_t primesInRange;
   uint32_t thisPrime, prevPrime;
   vector<uint64_t>  primes;
   vector<uint64_t>::iterator it;
   
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "primorial.pfgw";

   // Need to get ii_MinPrimorial and ii_MaxPrimorial from the input file
   if (is_InputTermsFileName.length() > 0)
      ProcessInputTermsFile(false);
   else 
   {
      if (ii_MaxPrimorial <= ii_MinPrimorial)
         FatalError("The value for -N must be greater than the value for -n");
   }
      
   // Count the list of primorial primes from FIRST_PRIMORIAL_PRIME to maxPrimorial
   primesInRange = primesieve::count_primes(FIRST_PRIMORIAL_PRIME + 1, ii_MaxPrimorial + 1);
   
   primes.clear();
   
   // Generate the list of primorial primes from FIRST_PRIMORIAL_PRIME to maxPrimorial
   primesieve::generate_n_primes(primesInRange, FIRST_PRIMORIAL_PRIME + 1, &primes);
   
   it = primes.begin();

   ip_PrimorialPrimes = (uint32_t *) xmalloc((primesInRange + 2) * sizeof(uint32_t));
   ip_PrimorialPrimeGaps = (uint16_t *) xmalloc((primesInRange + 2) * sizeof(uint16_t));
   
   ii_NumberOfPrimorialPrimes = 0;
   ii_BiggestGap = 0;
   prevPrime = FIRST_PRIMORIAL_PRIME;

   // Note that the max primorial is less than 2^32, so we can use a smaller datatype.
   while (it != primes.end())
   {
      thisPrime = (uint32_t) *it;
      it++;
      
      ip_PrimorialPrimes[ii_NumberOfPrimorialPrimes] = thisPrime;
      ip_PrimorialPrimeGaps[ii_NumberOfPrimorialPrimes] = (thisPrime - prevPrime);
      
      ii_NumberOfPrimorialPrimes++;

      if ((thisPrime - prevPrime) > ii_BiggestGap)
         ii_BiggestGap = (thisPrime - prevPrime);
         
      prevPrime = thisPrime;
   }   
   
   ip_PrimorialPrimes[ii_NumberOfPrimorialPrimes] = 0;
   ip_PrimorialPrimeGaps[ii_NumberOfPrimorialPrimes] = 0;

   iv_MinusTerms.resize(ii_MaxPrimorial - ii_MinPrimorial + 1);
   std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
   
   iv_PlusTerms.resize(ii_MaxPrimorial - ii_MinPrimorial + 1);
   std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
   
   if (is_InputTermsFileName.length() > 0)
   {
      il_TermCount = 0;
      ProcessInputTermsFile(true);
   }
   else
   {
      it = primes.begin();

      while (it != primes.end())
      {
         uint32_t thisPrime = (uint32_t) *it;
         it++;
         
         if (thisPrime >= ii_MinPrimorial && thisPrime <= ii_MaxPrimorial)
         {
            iv_PlusTerms[BIT(thisPrime)] = true;
            iv_MinusTerms[BIT(thisPrime)] = true;
            il_TermCount += 2;
         }
      }
   }

#ifdef HAVE_GPU_WORKERS
   SetMinGpuPrime(ii_MaxPrimorial + 1);
#endif

   FactorApp::ParentValidateOptions();
}

Worker *PrimorialApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      return new PrimorialGpuWorker(id, this);
#endif

   return new PrimorialWorker(id, this);
}

void PrimorialApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000], *pos;
   uint32_t primorial;
   int32_t  c;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
  if (memcmp(buffer, "ABC $a#$b", 9) && memcmp(buffer, "ABC $a#+$b", 10))
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   pos = strstr(buffer, "$b //");
   if (pos)
      if (sscanf(pos+16, "%" SCNu64"", &sieveLimit) == 1)
         SetMinPrime(sieveLimit);

   if (!haveBitMap)
      ii_MinPrimorial = ii_MaxPrimorial = 0;
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
   
      if (sscanf(buffer, "%u %d", &primorial, &c) != 2)
         FatalError("Line %s is malformed", buffer);

      if (!ii_MaxPrimorial)
         ii_MinPrimorial = ii_MaxPrimorial = primorial;
            
      if (haveBitMap)
      {
         if (c == -1)
         {
            iv_MinusTerms[primorial - ii_MinPrimorial] = true;
            il_TermCount++;
         }
         
         if (c == +1)
         {
            iv_PlusTerms[primorial - ii_MinPrimorial] = true;
            il_TermCount++;
         }
      }
      else
      {
         if (ii_MinPrimorial > primorial) ii_MinPrimorial = primorial;
         if (ii_MaxPrimorial < primorial) ii_MaxPrimorial = primorial;
      }
   }

   fclose(fPtr);
}

bool PrimorialApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t primorial;
   int32_t  c;
   
   if (sscanf(term, "%u#%d", &primorial, &c) != 2)
      FatalError("Could not parse term %s", term);
   
   if (primorial < ii_MinPrimorial || primorial > ii_MaxPrimorial)
      return false;
   
   if (!VerifyFactor(false, thePrime, primorial, c))
      return false;
   
   uint32_t bit = BIT(primorial);
   
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

void PrimorialApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint32_t termCount = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC $a#$b // Sieved to %" PRIu64"\n", largestPrime);

   for (uint32_t primorial=ii_MinPrimorial; primorial<=ii_MaxPrimorial; primorial++)
   {
      if (iv_MinusTerms[primorial - ii_MinPrimorial])
      {
         fprintf(termsFile, "%u -1\n", primorial);
         termCount++;
      }

      if (iv_PlusTerms[primorial - ii_MinPrimorial])
      {
         fprintf(termsFile, "%u +1\n", primorial);
         termCount++;
      }
   }

   fclose(termsFile);
   
   if (termCount != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termCount, il_TermCount);
   
   ip_FactorAppLock->Release();
}

void PrimorialApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%u <= primorial <= %u", ii_MinPrimorial, ii_MaxPrimorial);
}

bool PrimorialApp::ReportFactor(uint64_t primeFactor, uint32_t primorial, int32_t c, bool verifyFactor)
{
   uint32_t bit;
   bool     newFactor = false;
   
   if (primorial < ii_MinPrimorial || primorial > ii_MaxPrimorial)
      return false;

   if (verifyFactor)
      VerifyFactor(true, primeFactor, primorial, c);
   
   ip_FactorAppLock->Lock();

   bit = BIT(primorial);
   
   if (c == -1 && iv_MinusTerms[bit])
   {
      newFactor = true;
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(primeFactor, "%u#-1", primorial);
   }

   if (c == +1 && iv_PlusTerms[bit])
   {
      newFactor = true;
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(primeFactor, "%u#+1", primorial);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

bool  PrimorialApp::VerifyFactor(bool badFactorIsFatal, uint64_t primeFactor, uint32_t primorial, int32_t theC)
{
   uint64_t rem = FIRST_PRIMORIAL;
   uint32_t i = 0;
   
   fpu_push_1divp(primeFactor);
   
   for (i=0; ; i++)
   {
      rem = fpu_mulmod(rem, ip_PrimorialPrimes[i], primeFactor);

      if (ip_PrimorialPrimes[i] == primorial)
         break;
   }
  
   fpu_pop();

   if (theC == -1 && rem == +1)
      return true;
      
   if (theC == +1 && rem == primeFactor-1)
      return true;
       
   if (badFactorIsFatal)
      FatalError("%" PRIu64" is not a factor of %u#%+d (%llu)", primeFactor, primorial, theC, rem);

   return false;
}