/* SmarandacheApp.cpp -- (C) Mark Rodenkirch, January 2022

   This program finds factors of Smarandache numbers.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include "SmarandacheApp.h"
#include "SmarandacheWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../core/Parser.h"

#ifdef HAVE_GPU_WORKERS
#include "SmarandacheGpuWorker.h"
#define APP_NAME        "smsievecl"
#else
#define APP_NAME        "smsieve"
#endif

#define APP_VERSION     "1.0"

#define BIT(n)          ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the CPUSieve framework for other applications.
App *get_app(void)
{
   return new SmarandacheApp();
}

SmarandacheApp::SmarandacheApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of Smarandache numbers");
   SetLogFileName("smsieve.log");

   ii_MinN = 0;
   ii_MaxN = 0;
   ii_CpuWorkSize = 10000;
   
   // We'll remove all even terms manually
   SetAppMinPrime(99*99);

   SetAppMaxPrime(PMAX_MAX_52BIT);
   
#ifdef HAVE_GPU_WORKERS
   ii_MaxGpuSteps = 1000000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 10;
#endif
}

void SmarandacheApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=M           maximum n to search\n");

#ifdef HAVE_GPU_WORKERS
   printf("-S --step=S           max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  SmarandacheApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
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

parse_t SmarandacheApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'n':
         status = Parser::Parse(arg, 100000, 9999999, ii_MinN);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 100000, 9999999, ii_MaxN);
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

void SmarandacheApp::ValidateOptions(void)
{
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[50];
      
      sprintf(fileName, "smarandache.pfgw");
      
      is_OutputTermsFileName = fileName;
   }
   
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
      
      il_TermCount = (ii_MaxN - ii_MinN + 1);
   }

   SetMinGpuPrime(ii_MaxN + 1);

   FactorApp::ParentValidateOptions();

   // The testing routine is optimized to test 4 primes at a time.
   while (ii_CpuWorkSize % 4 > 0)
      ii_CpuWorkSize++;
}

Worker *SmarandacheApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      return new SmarandacheGpuWorker(id, this);
#endif
   
   return new SmarandacheWorker(id, this);
}

terms_t *SmarandacheApp::GetTerms(void)
{
   terms_t *terms = (terms_t *) xmalloc(sizeof(terms_t));
   uint64_t idx = 0;
   
   terms->termList = (uint32_t *) xmalloc(il_TermCount * sizeof(uint32_t));
   
   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {      
      if (iv_Terms[BIT(n)])
      {
         terms->termList[idx] = n;
         idx++;
      }
   }
   
   terms->termCount = idx;
   return terms;
}

void SmarandacheApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t n;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
  if (memcmp(buffer, "ABC Sm($a)", 10))
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC Sm($a) // Sieved to %" SCNu64"", &sieveLimit) == 1)
      SetMinPrime(sieveLimit);

   if (!haveBitMap)
      ii_MinN = ii_MaxN = 0;
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
   
      if (sscanf(buffer, "%u", &n) != 1)
         FatalError("Line %s is malformed", buffer);

      if (n % 6 != 1)
         FatalError("Term %u is not == 1 (mod 6)", n);

      if (n<100000 || n>9999999)
         FatalError("Term %u must be between 100000 and 9999999", n);
      
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

bool SmarandacheApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t n;
   
   if (sscanf(term, "Sm(%u)", &n) != 1)
      FatalError("Could not parse term %s", term);
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
     
   if (!VerifyFactor(false, thePrime, n))
      return false;
   
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

void SmarandacheApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint64_t termsCounted = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC Sm($a) // Sieved to %" PRIu64"\n", largestPrime);

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

void SmarandacheApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%u <= n <= %u", ii_MinN, ii_MaxN);
}

bool SmarandacheApp::ReportFactor(uint64_t p, uint32_t n)
{
   uint32_t bit;
   bool     newFactor = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   ip_FactorAppLock->Lock();
   
   bit = BIT(n);
   
   if (iv_Terms[bit])
   {
      newFactor = true;
      iv_Terms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(p, "Sm(%u)", n);
   }
   
   ip_FactorAppLock->Release();
   
   VerifyFactor(true, p, n);
   
   return newFactor;
}

bool  SmarandacheApp::VerifyFactor(bool badFactorIsFatal, uint64_t p, uint32_t n)
{
   uint64_t rem = 1, i;
   char     buffer[100];
   
   fpu_push_1divp(p);
   
   for (i=2; i<=n; i++)
   {
      if (i < 10)
      {
         rem = fpu_mulmod(rem, 10, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }

      if (i < 100)
      {
         rem = fpu_mulmod(rem, 100, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
      
      if (i < 1000)
      {
         rem = fpu_mulmod(rem, 1000, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
      
      if (i < 10000)
      {
         rem = fpu_mulmod(rem, 10000, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
      
      if (i < 100000)
      {
         rem = fpu_mulmod(rem, 100000, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
      
      if (i < 1000000)
      {
         rem = fpu_mulmod(rem, 1000000, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
      
      if (i < 10000000)
      {
         rem = fpu_mulmod(rem, 10000000, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
      
      if (i < 100000000)
      {
         rem = fpu_mulmod(rem, 100000000, p);
         rem += i;
         while (rem >= p)
            rem -= p;
         
         continue;
      }
   }
   
   fpu_pop();

   if (rem == 0)
      return true;

   sprintf(buffer, "%" PRIu64" is not a factor of not a factor of Sm(%u)", p, n);

   if (badFactorIsFatal)
      FatalError(buffer);
   else
      WriteToConsole(COT_OTHER, buffer);

   return false;
}