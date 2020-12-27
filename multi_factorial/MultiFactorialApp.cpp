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
#ifdef HAVE_GPU_WORKERS
#include "MultiFactorialGpuWorker.h"
#endif
#include "../core/Parser.h"

#define APP_NAME        "mfsieve"
#define APP_VERSION     "1.7"

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
   
   ip_FactorValidator = new MultiFactorialWorker(0, this);

   // We'll remove all even terms manually
   SetAppMinPrime(3);

   SetAppMaxPrime(PMAX_MAX_52BIT);
   
#ifdef HAVE_GPU_WORKERS
   ii_MaxGpuSteps = 100000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 10;
#endif
}

void MultiFactorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=M           maximum n to search\n");
   printf("-m --multifactorial=m multifactorial, e.g. x!m where m = 3 --> x!!! (default 1)\n");

#ifdef HAVE_GPU_WORKERS
   printf("-S --step=S           max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  MultiFactorialApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:m:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');
   AppendLongOpt(longOpts, "multifactorial", required_argument, 0, 'm');

#ifdef HAVE_GPU_WORKERS
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

void MultiFactorialApp::ValidateOptions(void)
{ 
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[50];
      
      if (ii_MultiFactorial == 1)
         sprintf(fileName, "factorial.pfgw");
      else
         sprintf(fileName, "mf_%d.pfgw", ii_MultiFactorial);
      
      is_OutputTermsFileName = fileName;
   }
   
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

   FactorApp::ParentValidateOptions();

   // The testing routine is optimized to test 4 primes at a time.
   while (ii_CpuWorkSize % 4 > 0)
      ii_CpuWorkSize++;
}

Worker *MultiFactorialApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      return new MultiFactorialGpuWorker(id, this);
#endif
   
   return new MultiFactorialWorker(id, this);
}

terms_t *MultiFactorialApp::GetTerms(void)
{
   terms_t *terms = (terms_t *) xmalloc(ii_MultiFactorial * sizeof(terms_t));
   double minPrime = (double) GetMinPrime();
   uint64_t bigN, bigMinN = ii_MinN;

   // Build an array by multiplying terms so that we can reduce the number of mulmods  in the main loop.
   // For example if multiFactorial = 1, minPrime = 1e6 and minN = 25:
   //   terms[0].termList[0] = 1 * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9
   //   terms[0].termList[1] = 10 * 11 * 12 * 13 * 14
   //   terms[0].termList[2] = 15 * 16 * 17 * 18
   //   terms[0].termList[3] = 19 * 20 * 21 * 22
   //   terms[0].termList[4] = 23 * 24
   // Note that each term < 1e6 (minPrime).
   // In short, instead of needing 24 mulmods to compute 25! (mod p) we will only need 5 mulmods.
   for (uint32_t mf=1; mf<=ii_MultiFactorial; mf++)
   {
      uint64_t idx = 0;
      
      terms[mf-1].mf = mf;
      terms[mf-1].termList = (uint64_t *) xmalloc(ii_MinN * sizeof(uint64_t)/2);
      
      terms[mf-1].termList[idx] = 1;
      for (uint32_t n=mf; n<ii_MinN; n+=mf)
      {
         if (terms[mf-1].termList[idx] * n > minPrime)
         {
            idx++;
            terms[mf-1].maxNForTerm = n;
            terms[mf-1].termList[idx] = n;
            continue;
         }

         terms[mf-1].maxNForTerm = n;
         terms[mf-1].termList[idx] *= n;
         
         bigN = n;
         if ((bigN * bigN) > bigMinN)
            break;
      }
   }
   
   return terms;
}

void MultiFactorialApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[200], *pos;
   uint32_t n;
   int32_t  c;
   uint64_t sieveLimit;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, 1000, fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
  if (memcmp(buffer, "ABC $a#!$b", 9) && memcmp(buffer, "ABC $a!+$b", 10) &&
      sscanf(buffer, "ABC $a!%d$b", &ii_MultiFactorial) != 1)
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   pos = strstr(buffer, "//");
   if (pos)
      if (sscanf(pos+13, "%" SCNu64"", &sieveLimit) == 1)
         SetMinPrime(sieveLimit);

   if (!haveBitMap)
      ii_MinN = ii_MaxN = 0;
   
   while (fgets(buffer, 1000, fPtr) != NULL)
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

bool MultiFactorialApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t n, mf;
   int32_t  c;
   
   if (sscanf(term, "%u!%u%d", &n, &mf, &c) != 3)
      FatalError("Could not parse term %s", term);

   if (mf != ii_MultiFactorial)
      FatalError("Expected multifactorial of %u in factor, but found %u", ii_MultiFactorial, mf);
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
  
   MultiFactorialWorker *mfWorker = (MultiFactorialWorker *) ip_FactorValidator;
   
   if (!mfWorker->VerifyFactor(false, thePrime, n, c))
   {
      if (mf == 1)
         WriteToConsole(COT_OTHER, "%" PRIu64" is not a factor of %u!%+d and was rejected", thePrime, n, c);
      else
         WriteToConsole(COT_OTHER, "%" PRIu64" is not a factor of %u!%u%+d and was rejected", thePrime, n, mf, c);
      
      return false;
   }
   
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
   uint32_t termCount = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC $a!%d$b // Sieved to %" PRIu64"\n", ii_MultiFactorial, largestPrime);

   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      if (iv_MinusTerms[n - ii_MinN])
      {
         fprintf(termsFile, "%d -1\n", n);
         termCount++;
      }

      if (iv_PlusTerms[n - ii_MinN])
      {
         fprintf(termsFile, "%d +1\n", n);
         termCount++;
      }
   }

   fclose(termsFile);
   
   if (termCount != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termCount, il_TermCount);
   
   ip_FactorAppLock->Release();
}

void MultiFactorialApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   if (ii_MultiFactorial > 1)
      sprintf(extraTtext, "%d <= n <= %d, multifactorial %u", ii_MinN, ii_MaxN, ii_MultiFactorial);
   else
      sprintf(extraTtext, "%d <= n <= %d, factorial", ii_MinN, ii_MaxN);
}

bool MultiFactorialApp::ReportFactor(uint64_t p, uint32_t n, int32_t c)
{
   uint32_t bit;
   bool     newFactor = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   ip_FactorAppLock->Lock();

   bit = BIT(n);
   
   if (c == -1 && iv_MinusTerms[bit])
   {
      newFactor = true;
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(p, "%u!%u-1", n, ii_MultiFactorial);
   }

   if (c == +1 && iv_PlusTerms[bit])
   {
      newFactor = true;
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(p, "%u!%u+1", n, ii_MultiFactorial);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

void MultiFactorialApp::ReportPrime(uint64_t p, uint32_t n, int32_t c)
{
      
   WriteToConsole(COT_OTHER, "%d!%d%+d is prime! (%" PRId64")", n, ii_MultiFactorial, c, p);

   WriteToLog("%d!%d%+d is prime! (%" PRId64")", n, ii_MultiFactorial, c, p);
}
