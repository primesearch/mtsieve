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
#include "../sieve/primesieve.hpp"

#define APP_NAME        "psieve"
#define APP_VERSION     "1.3"

#define BIT(n)          ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the CPUSieve framework for other applications.
App *get_app(void)
{
   return new PrimorialApp();
}

PrimorialApp::PrimorialApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of primorial");
   SetLogFileName("psieve.log");

   ii_MinN = 3;
   ii_MaxN = 0;
   ii_CpuWorkSize = 50000;
   
   // We'll remove all even terms manually
   SetAppMinPrime(3);
   
   // This is because the assembly code is using SSE to do the mulmods
   SetAppMaxPrime(PMAX_MAX_52BIT);
}

void PrimorialApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=M           maximum n to search\n");
}

void  PrimorialApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "n:N:";

   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');
}

parse_t PrimorialApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {         
      case 'n':
         status = Parser::Parse(arg, 3, 1000000000, ii_MinN);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 3, 1000000000, ii_MaxN);
         break;
   }

   return status;
}

void PrimorialApp::ValidateOptions(void)
{
   uint32_t primesInRange;
   vector<uint64_t>  primes;
   vector<uint64_t>::iterator it;
   
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "primorial.pfgw";
   
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
      if (ii_MaxN <= ii_MinN)
         FatalError("The value for -N must be greater than the value for -n");

      iv_MinusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
      
      iv_PlusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
      
      primesInRange = primesieve::count_primes(0, ii_MaxN + 1);
      
      primes.clear();
      
      // Generate primes for this worker
      primesieve::generate_n_primes(primesInRange, 0, &primes);
      
      it = primes.begin();
      
      while (it != primes.end())
      {
         if (*it >= ii_MinN && *it <= ii_MaxN)
         {
            iv_PlusTerms[BIT(*it)] = true;
            iv_MinusTerms[BIT(*it)] = true;
            il_TermCount += 2;
         }
         it++;
      }
   }

   FactorApp::ParentValidateOptions();
}

Worker *PrimorialApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{  
   return new PrimorialWorker(id, this);
}

void PrimorialApp::ProcessInputTermsFile(bool haveBitMap)
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
   
  if (memcmp(buffer, "ABC $a#$b", 9) && memcmp(buffer, "ABC $a#+$b", 10))
      FatalError("Line 1 is malformed in input file %s", is_InputTermsFileName.c_str());

   pos = strstr(buffer, "$b //");
   if (pos)
      if (sscanf(pos+16, "%" SCNu64"", &sieveLimit) == 1)
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

bool PrimorialApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t n;
   int32_t  c;
   
   if (sscanf(term, "%u#%d", &n, &c) != 2)
      FatalError("Could not parse term %s", term);
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
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

void PrimorialApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");
   uint32_t termCount = 0;

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC $a#$b // Sieved to %" PRIu64"\n", largestPrime);

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

void PrimorialApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%d <= primes <= %d", ii_MinN, ii_MaxN);
}

bool PrimorialApp::ReportFactor(uint64_t p, uint32_t n, int32_t c)
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
      
      LogFactor(p, "%u#-1", n);
   }

   if (c == +1 && iv_PlusTerms[bit])
   {
      newFactor = true;
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      
      LogFactor(p, "%u#+1", n);
   }
   
   ip_FactorAppLock->Release();
   
   return newFactor;
}

void PrimorialApp::ReportPrime(uint64_t p, uint32_t n, int32_t c)
{
   uint32_t bit;
   
   if (n < ii_MinN)
      return;
   
   if (n > ii_MaxN)
      return;

   ip_FactorAppLock->Lock();
   
   bit = BIT(n);
   
   if (c == -1 && iv_MinusTerms[bit])
   {	
      iv_MinusTerms[bit] = false;
      il_TermCount--;
   }

   if (c == +1 && iv_PlusTerms[bit])
   {	
      iv_PlusTerms[bit] = false;
      il_TermCount--;
   }
   
   WriteToConsole(COT_OTHER, "%d#%+d is prime! (%" PRId64")", n, c, p);

   WriteToLog("%d#%+d is prime! (%" PRId64")", n, c, p);
      
   ip_FactorAppLock->Release();
}
