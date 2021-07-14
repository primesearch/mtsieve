/* K1B2App.cpp -- (C) Mark Rodenkirch, November 2017

   Sieve for k*2^n+1 for a range of k and a range of n

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "../core/inline.h"
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "K1B2App.h"
#include "K1B2Worker.h"
#include "../x86_asm_ext/asm-ext-x86.h"

#define CMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

#define APP_NAME        "k1b2sieve"
#define APP_VERSION     "1.1"

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new K1B2App();
}

K1B2App::K1B2App(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of 2^n+c numbers for variable n and c");
   SetLogFileName("k1b2sieve.log");

   ii_MinN = 0;
   il_MaxC = 0;
   il_MaxC = 0;

   SetAppMinPrime(3);
}

void K1B2App::Help(void)
{
   FactorApp::ParentHelp();

   printf("-n --nmin=N           minimum n to search\n");
   printf("-N --nmax=N           maximum n to search\n");
   printf("-c --cmin=c           Minimum c to search\n");
   printf("-C --cmax=C           Maximum c to search\n");
}

void  K1B2App::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "c:C:n:N:";

   AppendLongOpt(longOpts, "nmin",              required_argument, 0, 'n');
   AppendLongOpt(longOpts, "nmax",              required_argument, 0, 'N');
   AppendLongOpt(longOpts, "cmin",              required_argument, 0, 'c');
   AppendLongOpt(longOpts, "cmax",              required_argument, 0, 'C');
}

parse_t K1B2App::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'c':
         status = Parser::Parse(arg, -CMAX_MAX, CMAX_MAX, il_MinC);
         break;

      case 'C':
         status = Parser::Parse(arg, -CMAX_MAX, CMAX_MAX, il_MaxC);
         break;
         
      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MinN);
         break;

      case 'N':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MaxN);
         break;
   }

   return status;
}

void K1B2App::ValidateOptions(void)
{
   uint32_t nCount;
   uint64_t cCount;
   
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "k1b2.pfgw";

   if (is_InputTermsFileName.length() > 0)
   {     
      ProcessInputTermsFile(false);
            
      nCount = ii_MaxN - ii_MinN + 1;
      cCount = il_MaxC - il_MinC + 1;
      
      iv_Terms.resize(nCount);
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         iv_Terms[n-ii_MinN].resize(cCount);
         std::fill(iv_Terms[n-ii_MinN].begin(), iv_Terms[n-ii_MinN].end(), false);
      }
            
      ProcessInputTermsFile(true);
   }
   else
   {
      if (il_MinC == 0)
         FatalError("kmin must be specified");

      if (il_MaxC == 0)
         FatalError("kmax must be specified");
      
      if (ii_MinN == 0)
         FatalError("nmin must be specified");

      if (ii_MaxN == 0)
         FatalError("nmax must be specified");
      
      if (il_MaxC <= il_MinC)
         FatalError("cmax must be greater than cmin");
      
      if (ii_MaxN < ii_MinN)
         FatalError("nmax must be greater than nmin");

      if (!(il_MinC & 1))
         il_MinC++;
      
      if (!(il_MaxC & 1))
         il_MaxC--;

      nCount = ii_MaxN - ii_MinN + 1;
      cCount = il_MaxC - il_MinC + 1;
      
      il_TermCount = 0;
      iv_Terms.resize(nCount);
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         iv_Terms[n-ii_MinN].resize(cCount);
         std::fill(iv_Terms[n-ii_MinN].begin(), iv_Terms[n-ii_MinN].end(), false);
         
         for (int64_t c=il_MinC; c<=il_MaxC; c+=2)
         {
            iv_Terms[n-ii_MinN][c-il_MinC] = true;
            il_TermCount++;
         }
      }
   }
   
   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(10000);
   
   SetMinGpuPrime(il_MaxC+1);
}

Worker *K1B2App::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      return new K1B2GpuWorker(id, this);
#endif

   return new K1B2Worker(id, this);
}

void K1B2App::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t n;
   int64_t  c;
   uint64_t lastPrime;
   bool     firstTerm = true;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
   if (!haveBitMap)
      il_MinC = il_MaxC = 0;
   
   if (!memcmp(buffer, "ABC ", 4))
   {
      if (sscanf(buffer, "ABC 2^$a$b // Sieved to %" SCNu64"", &lastPrime) != 1)
         FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());
      
      SetMinPrime(lastPrime);
   }
   else
      FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (sscanf(buffer, "%u %" SCNd64"", &n, &c) != 2)
         FatalError("Line %s is malformed", buffer);
    
      if (haveBitMap)
      {
         iv_Terms[n-ii_MinN][c-il_MinC] = true;
         il_TermCount++;
      }
      else
      {
         if (firstTerm)
         {
            ii_MinN = ii_MaxN = n;
            il_MinC = il_MaxC = c;
            firstTerm = false;
         }

         if (n < ii_MinN) ii_MinN = n;
         if (n > ii_MaxN) ii_MaxN = n;
         
         if (c < il_MinC) il_MinC = c;
         if (c > il_MaxC) il_MaxC = c;
      }
   }

   fclose(fPtr);
}

bool  K1B2App::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t n;
   int64_t  c;
   
   if (sscanf(term, "2^%u%" SCNd64"", &n, &c) != 2)
      FatalError("Could not parse term %s", term);

   if (n < ii_MinN || n > ii_MaxN)
      return false;
        
   if (c < il_MinC || c > il_MaxC)
      return false;
   
   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[n-ii_MinN][c-il_MinC])
   {
      iv_Terms[n-ii_MinN][c-il_MinC] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void K1B2App::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint32_t n;
   int64_t  c;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
   
   fprintf(termsFile, "ABC 2^$a$b // Sieved to %" PRIu64"\n", largestPrime);
   
   for (n=ii_MinN; n<=ii_MaxN; n++)
      for (c=il_MinC; c<=il_MaxC; c++)
      {
         if (iv_Terms[n-ii_MinN][c-il_MinC])
         {
            fprintf(termsFile, "%u %+" PRId64"\n", n, c);
            termsCounted++;
         }
      }

   fclose(termsFile);

   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void  K1B2App::GetExtraTextForSieveStartedMessage(char *extraText)
{ 
   sprintf(extraText, "%u <= n <= %u, %" PRId64" <= c <= %" PRId64", 2^n+c", ii_MinN, ii_MaxN, il_MinC, il_MaxC);
}

bool  K1B2App::ReportFactor(uint64_t p, uint32_t n, int64_t c)
{
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   if (c < il_MinC || c > il_MaxC)
      return false;
   
   bool removedTerm = false;
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();
      
   if (iv_Terms[n-ii_MinN][c-il_MinC])
   {
      iv_Terms[n-ii_MinN][c-il_MinC] = false;
      
      il_FactorCount++;
      il_TermCount--;
      
      LogFactor(p, "2^%u%+" PRId64"", n, c);
      removedTerm = true;
   }
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();

   return removedTerm;
}
