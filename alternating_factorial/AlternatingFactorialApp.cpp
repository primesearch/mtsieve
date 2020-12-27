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
#include "../core/Clock.h"
#include "AlternatingFactorialApp.h"
#include "AlternatingFactorialWorker.h"
#ifdef HAVE_GPU_WORKERS
#include "AlternatingFactorialGpuWorker.h"
#endif

#define APP_NAME        "afsieve"
#define APP_VERSION     "1.1"

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
   ip_FactorValidator = new AlternatingFactorialWorker(0, this);
   
   // Override the default
   ii_CpuWorkSize = 10000;

#ifdef HAVE_GPU_WORKERS
   ii_StepN = 10000;
#endif
}

void AlternatingFactorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=N           maximum n to search\n");

#ifdef HAVE_GPU_WORKERS
   printf("-s --step=s           n iterated per call to GPU\n");
#endif
}

void  AlternatingFactorialApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');

#ifdef HAVE_GPU_WORKERS
   shortOpts += "s:";
   
   AppendLongOpt(longOpts, "step_n",         required_argument, 0, 's');
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
         
#ifdef HAVE_GPU_WORKERS
      case 's':
         status = Parser::Parse(arg, 1, 1000000000, ii_StepN);
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

#ifdef HAVE_GPU_WORKERS
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
   char     buffer[200], *pos;
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

bool AlternatingFactorialApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t c;
      
   if (sscanf(term, "af(%u)", &c) != 1)
      FatalError("Could not parse term %s\n", term);

   if (c < ii_MinN || c > ii_MaxN)
      return false;
   
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

bool AlternatingFactorialApp::ReportFactor(int64_t p, uint32_t term)
{
   bool newFactor = false;
   
   if (term < ii_MinN || term > ii_MaxN)
      return false;
   
   ip_FactorAppLock->Lock();

   uint32_t bit = BIT(term);
   
   if (iv_Terms[bit])
   {
      newFactor = true;
      iv_Terms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(p, "af(%u)", term);
   }

   ip_FactorAppLock->Release();
   
   return newFactor;
}

void AlternatingFactorialApp::ReportPrime(int64_t p, uint32_t term)
{
   if (term < ii_MinN || term > ii_MaxN)
      return;

   ip_FactorAppLock->Lock();

   if (iv_Terms[BIT(term)])

   iv_Terms[BIT(term)] = false;

   WriteToConsole(COT_OTHER, "af(%d) is prime! (%" PRId64")", term, p);

   WriteToLog("af(%d) is prime! (%" PRId64")\n", term, p);

   il_TermCount--;

   ip_FactorAppLock->Release();
}

void AlternatingFactorialApp::WriteOutputTermsFile(uint64_t checkpointPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint32_t remaining = 0, bit;

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
   
   if (remaining != il_TermCount)
      FatalError("Terms expected != terms counted (%d != %d)", il_TermCount, remaining);
        
   ip_FactorAppLock->Release();
}
