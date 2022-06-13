/* SophieGermainApp.cpp -- (C) Mark Rodenkirch, July 2020

   Sieve for k*b^n-1 and k*b^(n+1)-1 for a range of k and a fixed b and n.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "../core/MpArith.h"
#include "SophieGermainApp.h"
#include "SophieGermainWorker.h"

#define APP_NAME        "sgsieve"
#define APP_VERSION     "1.3"

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
   SetBanner(APP_NAME " v" APP_VERSION ", a program to eliminate terms for Sophie-Germain prime searches for k*b^n-1 with variable k and fixed b and n");
   SetLogFileName("sgsieve.log");
   
   it_Format = FF_ABCD;
   
   il_MinK = 0;
   il_MaxK = 0;
   ii_Base = 0;
   ii_N = 0;
   ib_GeneralizedSearch = false;
 
   SetAppMinPrime(3);
   
   iv_Terms.clear();
}

void SophieGermainApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-b --base=b           Base to search\n");
   printf("-n --exp=n            Exponent to search\n");
   printf("-g --generalized      Multiply second term by b instead of by 2\n");
   printf("-f --format=f         Format of output file (D=ABCD (default), N=NEWPGEN)\n");
}

void  SophieGermainApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:b:n:gf:";

   AppendLongOpt(longOpts, "kmin",           required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",           required_argument, 0, 'K');
   AppendLongOpt(longOpts, "base",           required_argument, 0, 'b');
   AppendLongOpt(longOpts, "exp",            required_argument, 0, 'n');
   AppendLongOpt(longOpts, "generalized",    no_argument,       0, 'g');
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

      case 'b':
         status = Parser::Parse(arg, 2, NMAX_MAX, ii_Base);
         break;
         
      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_N);
         break;
      
      case 'g':
         ib_GeneralizedSearch = true;
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
      if (ib_GeneralizedSearch)
         WriteToConsole(COT_OTHER, "will not override -g, will be determined by input");
         
      ib_GeneralizedSearch = false;
      
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

      if (ii_Base == 0)
         FatalError("base must be specified");

      if (ii_N == 0)
         FatalError("exponent must be specified");

      if (ii_Base & 1)
      {
         if (il_MinK &1)
         {
            WriteToConsole(COT_OTHER, "kmin incremented by 1 as it must be even");
            il_MinK++;
         }

         if (il_MaxK &1)
         {
            WriteToConsole(COT_OTHER, "kmax decremented by 1 as it must be even");
            il_MaxK--;
         }
      }
      else
      {
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
      }
    
      il_TermCount = ((il_MaxK - il_MinK) >> 1) + 1;
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);      

   }
         
   if (is_OutputTermsFileName.length() == 0)
   {
      if (it_Format == FF_ABCD)
         is_OutputTermsFileName = "sg.abcd";
      else
      {
         if (!ib_GeneralizedSearch && ii_Base != 2)
            FatalError("newpgen format is not supported for this search");
         
         is_OutputTermsFileName = "sg.npg";
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
   uint32_t   b1, b2, n1, n2, n, m;
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
   
   if (sscanf(buffer, "ABCD $a*%u^%u-1 & %u*($a*%u^%u-1)+1 [%" SCNu64"] // Sieved to %" SCNu64"", &b1, &n1, &m, &b2, &n2, &k, &lastPrime) == 7) {
      if (b1 != b2)
         FatalError("Line 1 bases in ABCD input file %s do not match", is_InputTermsFileName.c_str());
      
      if (n1 != n2)
         FatalError("Line 1 exponents in ABCD input file %s do not match", is_InputTermsFileName.c_str());
   
      ii_Base = b1;
      ii_N = n1;
      format = FF_ABCD;
      if (m == 2)
         ib_GeneralizedSearch = false;
      else
      {
         if (m != b1)
            FatalError("Line 1 multiplier for the second term in ABCD input file %s is not supported", is_InputTermsFileName.c_str());

         ib_GeneralizedSearch = true;
      }

      if (haveBitMap)
      {
         bit = BIT(k);

         iv_Terms[bit] = true;
         il_TermCount++;
      }
      else 
         il_MinK = il_MaxK = k;
   }
   else if (sscanf(buffer, "%" SCNu64":S:0:%u:1034", &lastPrime, &ii_Base) == 2)
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

bool SophieGermainApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint64_t k;
   uint32_t b, n, m;
      
   if (sscanf(term, "%" SCNu64"*%u^%u-1", &k, &b, &n) != 2)
   {
      if (sscanf(term, "%u*(%" SCNu64"*%u^%u-1)+1", &m, &k, &b, &n) != 3)
         FatalError("Could not parse term %s", term);
      
      if (!ib_GeneralizedSearch && m != 2)
         FatalError("Expected multiplier of 2 in factor but found %u", b);
      
      if (ib_GeneralizedSearch && m != b)
         FatalError("Expected multiplier of %u in factor but found %u", ii_Base, b);
   }
   
   if (b != ii_Base)
      FatalError("Expected base %u in factor but found %u", ii_Base, b);
   
   if (n != ii_N && n != ii_N+1)
      FatalError("Expected n %u in factor but found %u", ii_N, n);
        
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
   uint64_t termsCounted = 0;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
      
   if (it_Format == FF_ABCD)
      termsCounted = WriteABCDTermsFile(largestPrime, termsFile);
   
   if (it_Format == FF_NEWPGEN)
      termsCounted = WriteNewPGenTermsFile(largestPrime, termsFile);
   
   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

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
   
   fprintf(termsFile, "ABCD $a*%u^%u-1 & %u*($a*%u^%u-1)+1 [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_Base, ii_N, (ib_GeneralizedSearch ? ii_Base : 2), ii_Base, ii_N, k, largestPrime);
   
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

   //   10 = k*2^n-1 & k*2^(n-1)+1 (per pfgw file newpgenformats.txt)
   // 1034 = k*2^n-1 & k*b^(n-1)+1 (per pfgw file newpgenformats.txt)
   fprintf(termsFile, "%" PRIu64":S:0:%u:%u\n", maxPrime, ii_Base, (ii_Base == 2 ? 10 : 1034));
      
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
   sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", first term t = k*%u^%u-1, next term is %u*t+1", il_MinK, il_MaxK, ii_Base, ii_N, (ib_GeneralizedSearch ? ii_Base : 2));
}

void  SophieGermainApp::ReportFactor(uint64_t theFactor, uint64_t k, bool firstOfPair, bool verifyFactor)
{
   if (ii_Base % theFactor == 0)
      return;
   
   if (verifyFactor)
      VerifyFactor(theFactor, k, firstOfPair);
   
   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();
      
   uint64_t bit = BIT(k);

   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      
      if (firstOfPair)
         LogFactor(theFactor, "%" PRIu64"*%u^%u-1", k, ii_Base, ii_N);
      else
         LogFactor(theFactor, "%u*(%" PRIu64"*%u^%u-1)+1", (ib_GeneralizedSearch ? ii_Base : 2), k, ii_Base, ii_N);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();
}

void  SophieGermainApp::VerifyFactor(uint64_t theFactor, uint64_t k, bool firstOfPair)
{
   MpArith  mp(theFactor);
   MpRes    pOne = mp.one();
   
   MpRes    resRem = mp.pow(mp.nToRes(ii_Base), ii_N);
   
   resRem = mp.mul(resRem, mp.nToRes(k));
   resRem = mp.sub(resRem, pOne);
   
   if (firstOfPair) 
   {
      if (resRem != mp.zero())
         FatalError("Invalid factor: %" PRIu64" is not a factor of not a factor of %" PRIu64"*%u^%u-1", theFactor, k, ii_Base, ii_N);
      
      return;
   }
   
   if (ib_GeneralizedSearch)
      resRem = mp.mul(resRem, mp.nToRes(ii_Base));
   else
      resRem = mp.mul(resRem, mp.nToRes(2));
      
   resRem = mp.add(resRem, pOne);
      
   if (resRem != mp.zero())
      FatalError("Invalid factor: %" PRIu64" is not a factor of not a factor of %u*(%" PRIu64"*%u^%u-1)+1", theFactor, (ib_GeneralizedSearch ? ii_Base : 2), k, ii_Base, ii_N);
   
}
