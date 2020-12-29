/* SophieGermainApp.cpp -- (C) Mark Rodenkirch, July 2020

   Sieve for k*b^n+1 and k*b^n-1 for a range of k and a fixed b and n.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "SophieGermainApp.h"
#include "SophieGermainWorker.h"

#define APP_NAME        "sgsieve"
#define APP_VERSION     "1.2"

#define NMAX_MAX        (1 << 31)

#define BIT(k)          (((k) - il_MinK) >> 1)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new SophieGermainApp();
}

SophieGermainApp::SophieGermainApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to eliminate terms for Sophie-Germain prime searches for base 2, fixed n and variable k");
   SetLogFileName("sgsieve.log");
   
   it_Format = FF_ABCD;
   
   il_MinK = 0;
   il_MaxK = 0;
   ii_N = 0;
 
   SetAppMinPrime(3);
   
   iv_Terms.clear();
}

void SophieGermainApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-n --exp=n            Exponent to search\n");
   printf("-f --format=f         Format of output file (D=ABCD (default), N=NEWPGEN)\n");
}

void  SophieGermainApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:n:f:";

   AppendLongOpt(longOpts, "kmin",           required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",           required_argument, 0, 'K');
   AppendLongOpt(longOpts, "exp",            required_argument, 0, 'n');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
}

parse_t SophieGermainApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'k':
         status = Parser::Parse(arg, 1, KMAX_MAX, il_MinK);
         break;

      case 'K':
         status = Parser::Parse(arg, 1, KMAX_MAX, il_MaxK);
         break;

      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_N);
         break;
         
      case 'f':
         char value;
         status = Parser::Parse(arg, "DN", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'D')
            it_Format = FF_ABCD;
         if (value == 'N')
            it_Format = FF_NEWPGEN;
         break;
   }

   return status;
}

void SophieGermainApp::ValidateOptions(void)
{
   if (it_Format == FF_UNKNOWN)
      FatalError("File format not valid, use D (ABCD) or N (NewPGen)");
   
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);

      iv_Terms.resize(((il_MaxK - il_MinK) >> 1) + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);

      ProcessInputTermsFile(true);
   }
   else
   {
      if (il_MinK == 0)
         FatalError("kmin must be specified");

      if (il_MaxK == 0)
         FatalError("kmax must be specified");
      
      if (il_MaxK <= il_MinK)
         FatalError("kmax must be greater than kmin");

      if (ii_N == 0)
         FatalError("exponent must be specified");

      if (!(il_MinK &1))
      {
         WriteToConsole(COT_OTHER, "kmin incremented by 1 as it must be odd");
         il_MinK++;
      }

      if (!(il_MaxK &1))
      {
         WriteToConsole(COT_OTHER, "kmax decremented by 1 as it must be odd");
         il_MaxK--;
      }
    
      il_TermCount = ((il_MaxK - il_MinK) >> 1) + 1;
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);      

   }
         
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "sg.abcd";

   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(10000);
}

Worker *SophieGermainApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that SophieGermain inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new SophieGermainWorker(id, this);

   return theWorker;
}

void SophieGermainApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE      *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char       buffer[1000];
   uint32_t   n1, n2, n;
   uint64_t   bit, k, diff, lastPrime;
   format_t   format = FF_UNKNOWN;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
   if (!haveBitMap)
   {
      ii_N = 0;
      il_MinK = il_MaxK = 0;
   }
   
   if (sscanf(buffer, "ABCD $a*2^%u+1 & 2*($a*2^%u+1)-1 [%" SCNu64"] // Sieved to %" SCNu64"", &n1, &n2, &k, &lastPrime) == 4) {
      if (n1 != n2)
         FatalError("Line 1 exponents in ABCD input file %s do not match", is_InputTermsFileName.c_str());
   
      ii_N = n1;
      format = FF_ABCD;

      if (haveBitMap)
      {
         bit = BIT(k);

         iv_Terms[bit] = true;
         il_TermCount++;
      }
      else 
         il_MinK = il_MaxK = k;
   }
   else if (sscanf(buffer, "%" SCNu64":C:0:2:5", &lastPrime) == 1)
   {
      format = FF_NEWPGEN;
   }
   else
      FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());     
   
   SetMinPrime(lastPrime);
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
      
      if (format == FF_ABCD)
      {
         if (sscanf(buffer, "%" SCNu64"", &diff) != 1)
            FatalError("Line %s is malformed", buffer);
         
         k += diff;

         if (haveBitMap)
         {
            bit = BIT(k);

            iv_Terms[bit] = true;
            il_TermCount++;
         }
         else
         {
            if (k == 0)
               FatalError("dd");
            
            if (il_MinK > k) il_MinK = k;
            if (il_MaxK < k) il_MaxK = k;
         }
      }
      else
      {
         if (sscanf(buffer, "%" SCNu64" %u", &k, &n) != 2)
            FatalError("Line %s is malformed", buffer);
         
         if (ii_N == 0)
         {
            ii_N = n;
            il_MinK = il_MaxK = k;
         }
         
         if (n != ii_N)
            FatalError("Only one n is supported");

         if (haveBitMap)
         {
            bit = BIT(k);

            iv_Terms[bit] = true;
            il_TermCount++;
         }
         else
         {
            if (il_MinK > k) il_MinK = k;
            if (il_MaxK < k) il_MaxK = k;
         }
      }
   }

   fclose(fPtr);
}

bool SophieGermainApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint64_t k;
   uint32_t n;
      
   if (sscanf(term, "%" SCNu64"*2^%u+1", &k, &n) != 2)
   {
      if (sscanf(term, "2*(%" SCNu64"*2^%u+1)-1", &k, &n) != 2)
         FatalError("Could not parse term %s", term);
   }
   
   if (n != ii_N && n != ii_N+1)
      FatalError("Expected n %u in factor but found %d", ii_N, n);
        
   if (k < il_MinK || k > il_MaxK)
      return false;
  
   uint64_t bit = k - il_MinK;
   
   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void SophieGermainApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t kCount = 0;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
      
   if (it_Format == FF_ABCD)
      kCount = WriteABCDTermsFile(largestPrime, termsFile);
   
   if (it_Format == FF_NEWPGEN)
      kCount = WriteNewPGenTermsFile(largestPrime, termsFile);
   
   fclose(termsFile);
   
   if (kCount != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", kCount, il_TermCount);
   
   ip_FactorAppLock->Release();
}

uint64_t SophieGermainApp::WriteABCDTermsFile(uint64_t largestPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0, previousK;
   uint64_t bit;
   
   for (k=il_MinK; k<=il_MaxK; k+=2)
   {
      bit = BIT(k);
   
      if (iv_Terms[bit])
         break;
   }
   
   if (k > il_MaxK)
      FatalError("No remaining terms");
   
   fprintf(termsFile, "ABCD $a*2^%d+1 & 2*($a*2^%d+1)-1 [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_N, ii_N, k, largestPrime);
   
   previousK = k;
   kCount = 1;
   k += 2;
   
   for (; k<=il_MaxK; k+=2)
   {
      bit = BIT(k);
      
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%" PRIu64"\n", k - previousK);
         previousK = k;
         kCount++;
      }
   }

   return kCount;
}

uint64_t SophieGermainApp::WriteNewPGenTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0;
   uint64_t bit;

   fprintf(termsFile, "%" PRIu64":C:0:2:5\n", maxPrime);
      
   k = il_MinK;
   
   for (; k<=il_MaxK; k+=2)
   {
      bit = BIT(k);
      
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%" PRIu64" %u\n", k, ii_N);
         kCount++;
      }
   }
   
   return kCount;
}

void  SophieGermainApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", k*2^%u+1 and 2*(k*2^%u+1)-1", il_MinK, il_MaxK, ii_N, ii_N);
}

bool  SophieGermainApp::ReportFactor(uint64_t p, uint64_t k, bool firstOfPair)
{
   bool     removedTerm = false;
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   uint64_t bit = BIT(k);

   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      removedTerm = true;
      
      if (firstOfPair)
         LogFactor(p, "%" PRIu64"*2^%u+1", k, ii_N);
      else
         LogFactor(p, "2*(%" PRIu64"*2^%u+1)-1", k, ii_N);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
   
   return removedTerm;
}
