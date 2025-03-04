/* SmarandacheWellinApp.cpp -- (C) Mark Rodenkirch, March 2025

   This program finds factors of SmarandacheWellin numbers.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include "SmarandacheWellinApp.h"
#include "SmarandacheWellinWorker.h"
#include "../core/MpArith.h"
#include "../core/Parser.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "SmarandacheWellinGpuWorker.h"
#define APP_NAME        "smwsievecl"
#else
#define APP_NAME        "smwsieve"
#endif

#define APP_VERSION     "1.0"

#define BIT(n)          ((n) - ii_MinN)

#define MAX_SUPPORTED_PRIME   1000000000

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the CPUSieve framework for other applications.
App *get_app(void)
{
   return new SmarandacheWellinApp();
}

SmarandacheWellinApp::SmarandacheWellinApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of Smarandache-Wellin numbers");
   SetLogFileName("smwsieve.log");

   ii_MinN = 0;
   ii_MaxN = 0;
   ii_CpuWorkSize = 10000;

   SetAppMinPrime(3);

   SetAppMaxPrime(PMAX_MAX_52BIT);
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuSteps = 1000000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}

void SmarandacheWellinApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=M           maximum n to search\n");

#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-S --step=S           max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  SmarandacheWellinApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
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

parse_t SmarandacheWellinApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'n':
         status = Parser::Parse(arg, 10000, MAX_SUPPORTED_PRIME, ii_MinN);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 10000, MAX_SUPPORTED_PRIME, ii_MaxN);
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

void SmarandacheWellinApp::ValidateOptions(void)
{
   std::vector<uint64_t>  primes;
   std::vector<uint64_t>::iterator it;


   if (is_InputTermsFileName.length() > 0)
      ProcessInputTermsFile(false);
   
   // Count the list of primes from 10 to maxN.  This is because we will start
   // with the first term > 2357, which is a Smarandache-Wellin prime
   uint32_t primesInRange = primesieve::count_primes(10, ii_MaxN + 1);
   
   primes.clear();
   
   // Generate the list of primes from 1 to maxN
   primesieve::generate_n_primes(primesInRange, 10, &primes);
   
   it = primes.begin();

   ip_Primes = (uint32_t *) xmalloc((primesInRange + 2), sizeof(uint32_t), "primes");
   ip_PrimeGaps = (uint16_t *) xmalloc((primesInRange + 2), sizeof(uint16_t), "primeGaps");
   
   ii_NumberOfPrimes = 0;

   uint64_t prevPrime = 7;
   ii_BiggestGap = 0;

   // Note that the max primorial is less than 2^32, so we can use a smaller datatype.
   while (it != primes.end())
   {
      uint64_t thisPrime = (uint32_t) *it;
      it++;
      
      if (thisPrime > ii_MaxN)
         break;
      
      ip_Primes[ii_NumberOfPrimes] = thisPrime;
      ip_PrimeGaps[ii_NumberOfPrimes] = (thisPrime - prevPrime);
      
      if ((thisPrime - prevPrime) > ii_BiggestGap)
         ii_BiggestGap = (thisPrime - prevPrime);
      
      prevPrime = thisPrime;
      
      ii_NumberOfPrimes++;
   }

   ip_Primes[ii_NumberOfPrimes] = 0;
   
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[50];
      
      snprintf(fileName, sizeof(fileName), "smarandache-wellin.pfgw");
      
      is_OutputTermsFileName = fileName;
   }
   
   if (is_InputTermsFileName.length() > 0)
   {   
      iv_Terms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);
      
      il_TermCount = 0;
      ProcessInputTermsFile(true);
   }
   else
   {
      if (ii_MaxN <= ii_MinN)
         FatalError("The value for -N must be greater than the value for -n");

      il_TermCount = 0;
      iv_Terms.resize(ii_MaxN - ii_MinN + 1);
      
      for (uint32_t i=0; i<ii_NumberOfPrimes; i++)
      {
         uint32_t n = ip_Primes[i];
         
         if (n >= ii_MinN && n <= ii_MaxN)
         {
            iv_Terms[BIT(n)] = true;
            il_TermCount++;
         }
      }
   }

   SetMinGpuPrime(10000);

   FactorApp::ParentValidateOptions();

   // The testing routine is optimized to test 4 primes at a time.
   while (ii_CpuWorkSize % 4 > 0)
      ii_CpuWorkSize++;
}

Worker *SmarandacheWellinApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#if defined(USE_OPENCL) || defined(USE_METAL)
   if (gpuWorker)
      return new SmarandacheWellinGpuWorker(id, this);
#endif
   
   return new SmarandacheWellinWorker(id, this);
}

void SmarandacheWellinApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t n;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
  if (memcmp(buffer, "ABC SmW($a)", 10))
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC SmW($a) // Sieved to %" SCNu64"", &sieveLimit) == 1)
      SetMinPrime(sieveLimit);

   if (!haveBitMap)
      ii_MinN = ii_MaxN = 0;
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
   
      if (sscanf(buffer, "%u", &n) != 1)
         FatalError("Line %s is malformed", buffer);

      if (n > MAX_SUPPORTED_PRIME)
         FatalError("Term %u cannot exceed %u", n, MAX_SUPPORTED_PRIME);
      
      if (ii_MaxN == 0)
         ii_MinN = ii_MaxN = n;
            
      if (haveBitMap)
      {
         iv_Terms[BIT(n)] = true;
         il_TermCount++;
      }
      else
      {
         if (ii_MinN > n) ii_MinN = n;
         if (ii_MaxN < n) ii_MaxN = n;
      }
   }

   fclose(fPtr);
}

void  SmarandacheWellinApp::FillTerms(uint8_t *terms)
{
   memset(terms, 0, ii_NumberOfPrimes);
   uint32_t term;

   for (uint32_t idx=0; idx<ii_NumberOfPrimes; idx++)
   {
      term = ip_Primes[idx];
      
      if (term >= ii_MinN)
         terms[idx] = (iv_Terms[BIT(term)] ? 1 : 0);
   }
}

bool SmarandacheWellinApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t n;
   
   if (sscanf(term, "SmW(%u)", &n) != 1)
      FatalError("Could not parse term %s", term);
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
     
   VerifyFactor(thePrime, n);
   
   uint64_t bit = BIT(n);
   
   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[bit])
   {	
      iv_Terms[bit] = false;
      il_TermCount--;
      return true;
   }
      
   return false;
}

void SmarandacheWellinApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint64_t termsCounted = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC SmW($a) // Sieved to %" PRIu64"\n", largestPrime);

   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      if (iv_Terms[BIT(n)])
      {
         fprintf(termsFile, "%d\n", n);
         termsCounted++;
      }
   }

   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void SmarandacheWellinApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   snprintf(extraText, maxTextLength, "%u <= n <= %u", ii_MinN, ii_MaxN);
}

bool SmarandacheWellinApp::ReportFactor(uint64_t theFactor, uint32_t n)
{
   uint32_t bit;
   bool     newFactor = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
      
   VerifyFactor(theFactor, n);
   
   ip_FactorAppLock->Lock();
   
   bit = BIT(n);
   
   if (iv_Terms[bit])
   {
      newFactor = true;
      iv_Terms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(theFactor, "SmW(%u)", n);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

void  SmarandacheWellinApp::VerifyFactor(uint64_t theFactor, uint32_t n)
{   
   MpArith mp(theFactor);
   MpRes   res = mp.nToRes(2357);
   MpRes   mp1e2 = mp.nToRes(100);
   MpRes   mp1e3 = mp.nToRes(1000);
   MpRes   mp1e4 = mp.nToRes(10000);
   MpRes   mp1e5 = mp.nToRes(100000);
   MpRes   mp1e6 = mp.nToRes(1000000);
   MpRes   mp1e7 = mp.nToRes(10000000);
   MpRes   mp1e8 = mp.nToRes(100000000);
   MpRes   mp1e9 = mp.nToRes(1000000000);

   for (uint32_t i=0; i<ii_NumberOfPrimes; i++)
   {
      if (ip_Primes[i] < 100)
         res = mp.mul(res, mp1e2);
      else if (ip_Primes[i] < 1000)
         res = mp.mul(res, mp1e3);
      else if (ip_Primes[i] < 10000)
         res = mp.mul(res, mp1e4);
      else if (ip_Primes[i] < 100000)
         res = mp.mul(res, mp1e5);
      else if (ip_Primes[i] < 1000000)
         res = mp.mul(res, mp1e6);
      else if (ip_Primes[i] < 10000000)
         res = mp.mul(res, mp1e7);
      else if (ip_Primes[i] < 100000000)
         res = mp.mul(res, mp1e8);
      else if (ip_Primes[i] < 1000000000)
         res = mp.mul(res, mp1e9);
      else
         FatalError("How did this happen?  n is too large");
      
      res = mp.add(res, mp.nToRes(ip_Primes[i]));
      
      if (ip_Primes[i] == n)
         break;
   }
   
   if (res == mp.zero())
      return;

   FatalError("%" PRIu64" is not a factor of not a factor of SmW(%u)", theFactor, n);
}