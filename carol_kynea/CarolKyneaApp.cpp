/* CarolKyneaApp.cpp -- (C) Mark Rodenkirch, July 2017

   CarolKynea application

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <time.h>
#include <cinttypes>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "../core/MpArith.h"
#include "CarolKyneaApp.h"
#include "CarolKyneaWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "CarolKyneaGpuWorker.h"
#define APP_NAME        "cksievecl"
#else
#define APP_NAME        "cksieve"
#endif

#define APP_VERSION     "1.4"

#define BIT(n)          ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new CarolKyneaApp();
}

CarolKyneaApp::CarolKyneaApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of (b^n+/-1)^2-2 numbers");
   SetLogFileName("cksieve.log");

   ii_Base = 0;
   ii_MinN = 1;
   ii_MaxN = 0;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}

void CarolKyneaApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-b --base=b           base to search\n");
   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=M           maximum n to search\n");
#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  CarolKyneaApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "b:n:N:";

   AppendLongOpt(longOpts, "minb",           required_argument, 0, 'b');
   AppendLongOpt(longOpts, "minn",           required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",           required_argument, 0, 'N');

#if defined(USE_OPENCL) || defined(USE_METAL)
   shortOpts += "M:";
   
   AppendLongOpt(longOpts, "maxfactors",     required_argument, 0, 'M');
#endif
}

parse_t CarolKyneaApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'b':
         status = Parser::Parse(arg, 2, 1000000000, ii_Base);
         break;
         
      case 'n':
         status = Parser::Parse(arg, 1, 1000000000, ii_MinN);
         break;
         
      case 'N':
         status = Parser::Parse(arg, 2, 1000000000, ii_MaxN);
         break;
         
#if defined(USE_OPENCL) || defined(USE_METAL)
      case 'M':
         status = Parser::Parse(arg, 10, 1000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void CarolKyneaApp::ValidateOptions(void)
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
      if (ii_Base == 0)
         FatalError("Base is required");
         
      if (ii_MaxN <= ii_MinN)
         FatalError("The value for -N must be greater than the value for -n");

      iv_MinusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), true);
      
      iv_PlusTerms.resize(ii_MaxN - ii_MinN + 1);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), true);
      
      il_TermCount = 2 * (ii_MaxN - ii_MinN + 1);

      // (2^1+1)-2 =  1, so remove that term
      if (ii_MinN == 1 && ii_Base == 2)
      {
         iv_PlusTerms[BIT(1)] = false;
         il_TermCount--;
      }
      
      // (2^1-1)-2 = -1, so remove that term
      // (3^1-1)-2 = -0, so remove that term
      // (4^1-1)-2 =  1, so remove that term
      if (ii_MinN == 1 && ii_Base <= 4)
      {
         iv_MinusTerms[BIT(1)] = false;
         il_TermCount--;
      }
   }

   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[50];
      
      sprintf(fileName, "ck_%u.pfgw", ii_Base);
      
      is_OutputTermsFileName = fileName;
   }

   FactorApp::ParentValidateOptions();

   SetMaxPrimeForSingleWorker(100000);
   
   // The GPU code for this sieve will not support primes lower than this.
   SetMinGpuPrime(10000);
}

Worker *CarolKyneaApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#if defined(USE_OPENCL) || defined(USE_METAL)  
   if (gpuWorker)
      theWorker = new CarolKyneaGpuWorker(id, this);
   else
#endif
      theWorker = new CarolKyneaWorker(id, this);

   return theWorker;
}

void CarolKyneaApp::ProcessInputTermsFile(bool haveBitMap)
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
   
   if (sscanf(buffer, "ABC (%u^$a$b)^2-2", &ii_Base) != 1)
      FatalError("Line 1 is not a valid ABC line in input file %s", is_InputTermsFileName.c_str());

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
   
      if (sscanf(buffer, "%u %d", &n, &c) != 2)
         FatalError("Line %s is malformed", buffer);

      if (!ii_MaxN)
         ii_MinN = ii_MaxN = n;

      if (haveBitMap)
      {
         if (c == +1)
         {
            iv_PlusTerms[BIT(n)] = true;
            il_TermCount++;
         }
         
         if (c == -1)
         {
            iv_MinusTerms[BIT(n)] = true;
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

bool CarolKyneaApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint32_t b, n;
   int32_t  c;
      
   if (sscanf(term, "(%u^%u%d)^2-2", &b, &n, &c) != 3)
      FatalError("Could not parse term %s", term);

   if (b != ii_Base)
      FatalError("Expected base %u in factor but found base %u", ii_Base, b);
     
   if (c != +1 && c != -1)
      FatalError("Expected c to be +1 or -1 in factor but found %d", c);
        
   if (n < ii_MinN || n > ii_MaxN)
      return false;

   uint32_t bit = BIT(n);
   
   // No locking is needed because the Workers aren't running yet
   if (c == +1 && iv_MinusTerms[bit])
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

void CarolKyneaApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint32_t n;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC (%u^$a$b)^2-2 // Sieved to %" PRIu64"\n", ii_Base, largestPrime);

   for (n=ii_MinN; n<=ii_MaxN; n++)
   {         
      if (iv_PlusTerms[n-ii_MinN])
      {
         fprintf(termsFile, "%u +1\n", n);
         termsCounted++;
      }
      
      if (iv_MinusTerms[n-ii_MinN])
      {
         fprintf(termsFile, "%u -1\n", n);
         termsCounted++;
      }
   }

   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void CarolKyneaApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%u <= n <= %u, (%u^n+/-1)^2-2", ii_MinN, ii_MaxN, ii_Base);
}

bool CarolKyneaApp::ReportFactor(uint64_t theFactor, uint32_t n, int32_t c)
{
   uint32_t bit;
   bool     removedTerm = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
      
   ip_FactorAppLock->Lock();

   bit = BIT(n);
   
   if (c == -1 && iv_MinusTerms[bit])
   {	
      iv_MinusTerms[bit] = false;
      removedTerm = true;
      
      il_TermCount--;
      il_FactorCount++;
      
      if (IsPrime(theFactor, n, c))
      {
         WriteToConsole(COT_OTHER, "(%u^%u-1)^2-2 is prime! (%" PRIu64")", ii_Base, n, theFactor);

         WriteToLog("(%u^%u-1)^2-2 is prime! (%" PRIu64")", ii_Base, n, theFactor);
      }
      else
         LogFactor(theFactor, "(%u^%u-1)^2-2", ii_Base, n);
      
      VerifyFactor(theFactor, n, c);
   }

   if (c == +1 && iv_PlusTerms[bit])
   {	
      iv_PlusTerms[bit] = false;
      removedTerm = true;
      
      il_TermCount--;
      il_FactorCount++;
      
      
      if (IsPrime(theFactor, n, c))
      {
         WriteToConsole(COT_OTHER, "(%u^%u+1)^2-2 is prime! (%" PRIu64")", ii_Base, n, theFactor);

         WriteToLog("(%u^%u+1)^2-2 is prime! (%" PRIu64")", ii_Base, n, theFactor);
      }
      else
         LogFactor(theFactor, "(%u^%u+1)^2-2", ii_Base, n);
      
      VerifyFactor(theFactor, n, c);
   }
   
   ip_FactorAppLock->Release();
   
   return removedTerm;   
}

bool CarolKyneaApp::IsPrime(uint64_t theFactor, uint32_t n, int32_t c)
{
   if (n >= 62)
      return false;
   
   uint64_t term = 1;
   uint64_t maxBeforeOverflow = (1 << 31);
   uint32_t i = 0;
   bool isFactor = false;
   
   do
   {
      term *= ii_Base;
      i++;
      
      if (term + c > maxBeforeOverflow)
         isFactor = true;
      
      if (((term + c) * (term + c)) > theFactor + 2)
         isFactor = true;
   } while (i < n && !isFactor);

   return !isFactor;
}

void CarolKyneaApp::VerifyFactor(uint64_t theFactor, uint32_t n, int32_t c)
{   
   if (n < ii_MinN || n > ii_MaxN)
      return;
   
   MpArith  mp(theFactor);
   MpRes    res = mp.pow(mp.nToRes(ii_Base), n);

   if (c == 1)
      res = mp.add(res, mp.one());
   else
      res = mp.sub(res, mp.one());
   
   res = mp.mul(res, res);
   
   uint64_t rem = mp.resToN(res);
   
   if (rem == 2)
      return;
      
   FatalError("%" PRIu64" is not a factor of (%u^%u%+d)^2-2", theFactor, ii_Base, n, c);
}
