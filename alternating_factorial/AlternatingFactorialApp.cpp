/* AlternatingFactorialApp.cpp -- (C) Mark Rodenkirch, July 2017

   AlternatingFactorial application

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include <cinttypes>
#include "../core/Parser.h"

#include "../core/MpArith.h"

#include "AlternatingFactorialApp.h"
#include "AlternatingFactorialWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "AlternatingFactorialGpuWorker.h"
#endif

#define APP_NAME        "afsieve"
#define APP_VERSION     "1.2"

#define BIT(n)          ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new AlternatingFactorialApp();
}

AlternatingFactorialApp::AlternatingFactorialApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of alternating factorials");
   SetLogFileName("afsieve.log");

   ii_MinN = 2;
   ii_MaxN = 0;
   
   // This is because the assembly code is using SSE to do the mulmods
   SetAppMaxPrime(PMAX_MAX_52BIT);
   
   // Override the default
   ii_CpuWorkSize = 10000;

#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuSteps = 100000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 10;
#endif
}

void AlternatingFactorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=N           maximum n to search\n");

#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-S --step=S           max steps iterated per call to GPU (default %d)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  AlternatingFactorialApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');

#if defined(USE_OPENCL) || defined(USE_METAL)
   shortOpts += "S:M:";
   
   AppendLongOpt(longOpts, "steps",             required_argument, 0, 'S');
   AppendLongOpt(longOpts, "maxfactors",        required_argument, 0, 'M');
#endif
}

parse_t AlternatingFactorialApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'n':
         status = Parser::Parse(arg, 2, 1000000000, ii_MinN);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 2, 1000000000, ii_MaxN);
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

void AlternatingFactorialApp::ValidateOptions(void)
{
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "af.terms";

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      iv_Terms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);

      il_TermCount = 0;
      ProcessInputTermsFile(true);
   }
   else
   {     
      if (ii_MaxN <= ii_MinN)
         FatalError("The value for -N must be greater than the value for -n");
      
      iv_Terms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);
      
      il_TermCount = ii_MaxN - ii_MinN + 1;
   }

   FactorApp::ParentValidateOptions();

   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
}

Worker *AlternatingFactorialApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#if defined(USE_OPENCL) || defined(USE_METAL)
   if (gpuWorker)
      theWorker = new AlternatingFactorialGpuWorker(id, this);
   else
#endif
      theWorker = new AlternatingFactorialWorker(id, this);
   
   return theWorker;
}

void AlternatingFactorialApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000], *pos;
   uint32_t n;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
  if (memcmp(buffer, "ABC af($a)", 10))
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
   
      if (sscanf(buffer, "%d", &n) != 1)
         FatalError("Line %s is malformed", buffer);

      if (!ii_MaxN)
         ii_MinN = ii_MaxN = n;

      if (haveBitMap)
      {
         il_TermCount++;
         iv_Terms[BIT(n)] = true;
      }
      else
      {
         if (ii_MinN > n) ii_MinN = n;
         if (ii_MaxN < n) ii_MaxN = n;
      }
   }

   fclose(fPtr);
}

bool AlternatingFactorialApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint32_t c;
      
   if (sscanf(term, "af(%u)", &c) != 1)
      FatalError("Could not parse term %s\n", term);

   if (c < ii_MinN || c > ii_MaxN)
      return false;
   
   VerifyFactor(theFactor, c);
   
   uint64_t bit = BIT(c);
   
   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void AlternatingFactorialApp::GetExtraTextForSieveStartedMessage(char *extraText)
{
   sprintf(extraText, "%d <= n <= %d", ii_MinN, ii_MaxN);
}

bool AlternatingFactorialApp::ReportFactor(uint64_t theFactor, uint32_t term)
{
   bool newFactor = false;
   
   if (term < ii_MinN || term > ii_MaxN)
      return false;
   
   VerifyFactor(theFactor, term);
   
   ip_FactorAppLock->Lock();

   uint32_t bit = BIT(term);
   
   if (iv_Terms[bit])
   {
      newFactor = true;
      iv_Terms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(theFactor, "af(%u)", term);
   }

   ip_FactorAppLock->Release();
   
   return newFactor;
}

void AlternatingFactorialApp::WriteOutputTermsFile(uint64_t checkpointPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint32_t remaining = 0, bit;
   uint32_t termCount = (uint32_t) il_TermCount;

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();

   fprintf(termsFile, "ABC af($a) // Sieved to %" PRIu64"\n", checkpointPrime);

   bit = BIT(ii_MinN);
   
   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%d\n", n);
         remaining++;
      }
      
      bit++;
   }

   fclose(termsFile);
   
   if (remaining != termCount)
      FatalError("Terms expected != terms counted (%u != %u)", termCount, remaining);
        
   ip_FactorAppLock->Release();
}

void  AlternatingFactorialApp::VerifyFactor(uint64_t theFactor, uint32_t term)
{
   MpArith  mp(theFactor);
   MpRes    mpN = mp.one();
   MpRes    mpFn = mp.one();
   MpRes    mpAfn = mp.one();
   
   for (uint32_t n=1; n<term; n++)
   {
      mpN = mp.add(mpN, mp.one());
      
      // Compute n!
      mpFn = mp.mul(mpFn, mpN);

      // Compute n! - af(n-1)!      
      mpAfn = mp.sub(mpFn, mpAfn);
   }

   uint64_t rem = mp.resToN(mpAfn);
   
   bool isValid = (rem == 0);
   
   if (isValid)
      return;
   
   FatalError("%" PRIu64" is not a factor of af(%u)", theFactor, term);
}
