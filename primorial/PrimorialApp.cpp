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

#include "../core/Parser.h"
#include "../core/MpArith.h"
#include "../sieve/primesieve.hpp"

#if defined(USE_OPENCL)
#include "PrimorialGpuWorker.h"
#define APP_NAME        "psievecl"
#elif defined(USE_METAL)
#include "PrimorialGpuWorker.h"
#define APP_NAME        "psievemtl"
#else
#define APP_NAME        "psieve"
#endif

// Since we have no explicit checks for p(n)+1 or p(n)-1 to be prime
// we want to ensure that p(n) > 124 bits.  All primorials below this are known
// so there isn't any value in sieving that range.
#define MIN_PRIMORIAL   100
#define MAX_PRIMORIAL   1000000000

#define APP_VERSION     "1.6"

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
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuSteps = 50000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}

void PrimorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum primorial to search\n");
   printf("-N --maxn=M           maximum primorial to search\n");

#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-S --step=S           max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  PrimorialApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');

#if defined(USE_OPENCL) || defined(USE_METAL)
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
         status = Parser::Parse(arg, MIN_PRIMORIAL, MAX_PRIMORIAL, ii_MinPrimorial);
         break;
         
      case 'N':
         status = Parser::Parse(arg, MIN_PRIMORIAL, MAX_PRIMORIAL, ii_MaxPrimorial);
         break;

#if defined(USE_OPENCL) || defined(USE_METAL)
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
   std::vector<uint64_t>  primes;
   std::vector<uint64_t>::iterator it;
   
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
      
   // Count the list of primorial primes from FIRST_PRIMORIAL_PRIME+1  to maxPrimorial
   primesInRange = primesieve::count_primes(FIRST_PRIMORIAL_PRIME + 1, ii_MaxPrimorial + 1);
   
   primes.clear();
   
   // Generate the list of primorial primes from FIRST_PRIMORIAL_PRIME+1 to maxPrimorial
   primesieve::generate_n_primes(primesInRange, FIRST_PRIMORIAL_PRIME + 1, &primes);
   
   it = primes.begin();

   ip_PrimorialPrimes = (uint32_t *) xmalloc((primesInRange + 2), sizeof(uint32_t), "primes");
   ip_PrimorialPrimeGaps = (uint16_t *) xmalloc((primesInRange + 2), sizeof(uint16_t), "primeGaps");
   
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

#if defined(USE_OPENCL) || defined(USE_METAL)
   SetMinGpuPrime(ii_MaxPrimorial + 1);
#endif

   FactorApp::ParentValidateOptions();
}

Worker *PrimorialApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#if defined(USE_OPENCL) || defined(USE_METAL)
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
   bool     skipped = false;

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

      if (primorial < MIN_PRIMORIAL || primorial > MAX_PRIMORIAL)
      {
         skipped = true;
         continue;
      }

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

   if (skipped && haveBitMap)
      WriteToConsole(COT_OTHER, "Igmoring primorials < %u or > %u as all primorials below this are known", MIN_PRIMORIAL, MAX_PRIMORIAL);
}

bool PrimorialApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t primorial;
   int32_t  c;
   
   if (sscanf(term, "%u#%d", &primorial, &c) != 2)
      FatalError("Could not parse term %s", term);
   
   if (primorial < ii_MinPrimorial || primorial > ii_MaxPrimorial)
      return false;
   
   VerifyFactor(thePrime, primorial, c);
   
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
   uint64_t termsCounted = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC $a#$b // Sieved to %" PRIu64"\n", largestPrime);

   for (uint32_t primorial=ii_MinPrimorial; primorial<=ii_MaxPrimorial; primorial++)
   {
      if (iv_MinusTerms[primorial - ii_MinPrimorial])
      {
         fprintf(termsFile, "%u -1\n", primorial);
         termsCounted++;
      }

      if (iv_PlusTerms[primorial - ii_MinPrimorial])
      {
         fprintf(termsFile, "%u +1\n", primorial);
         termsCounted++;
      }
   }

   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);
   
   ip_FactorAppLock->Release();
}

void PrimorialApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   snprintf(extraText, maxTextLength, "%u <= primorial <= %u", ii_MinPrimorial, ii_MaxPrimorial);
}

bool PrimorialApp::ReportFactor(uint64_t theFactor, uint32_t primorial, int32_t c)
{
   uint32_t bit;
   bool     newFactor = false;
   
   if (primorial < ii_MinPrimorial || primorial > ii_MaxPrimorial)
      return false;

   VerifyFactor(theFactor, primorial, c);
   
   ip_FactorAppLock->Lock();

   bit = BIT(primorial);
   
   if (c == -1 && iv_MinusTerms[bit])
   {
      newFactor = true;
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(theFactor, "%u#-1", primorial);
   }

   if (c == +1 && iv_PlusTerms[bit])
   {
      newFactor = true;
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(theFactor, "%u#+1", primorial);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

void  PrimorialApp::VerifyFactor(uint64_t theFactor, uint32_t primorial, int32_t c)
{
   MpArith  mp(theFactor);
   MpRes    pOne = mp.one();
   MpRes    mOne = mp.sub(mp.zero(), pOne);
   MpRes    res = mp.nToRes(FIRST_PRIMORIAL);

   for (uint32_t i=0; ; i++)
   {
      res = mp.mul(res, mp.nToRes(ip_PrimorialPrimes[i]));

      if (ip_PrimorialPrimes[i] == primorial)
         break;
   }
   
   if (c == -1 && res != pOne)
      FatalError("%" PRIu64" is not a factor of %u#-1", theFactor, primorial);
      
   if (c == +1 && res != mOne)
      FatalError("%" PRIu64" is not a factor of %u#+1", theFactor, primorial);
}