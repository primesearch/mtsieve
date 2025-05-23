/* FixedKBNApp.cpp -- (C) Mark Rodenkirch, November 2017

   Sieve for k*2^n+1 for a range of k and a range of n

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "FixedKBNApp.h"
#include "FixedKBNWorker.h"
#include "../core/MpArith.h"

#define APP_NAME        "fkbnsieve"
#define APP_VERSION     "1.6.4"

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new FixedKBNApp();
}

FixedKBNApp::FixedKBNApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of k*b^n+c numbers for fixed k, b, and n and variable c");
   SetLogFileName("fkbnsieve.log");

   SetAppMinPrime(3);
   is_Sequence = "";

   il_MinC = 0;
   il_MaxC = 0;
}

void FixedKBNApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-c --cmin=c           Minimum c to search\n");
   printf("-C --cmax=C           Maximum c to search\n");
   printf("-s --sequence=s       Sequence to find factors of in form k*b^n+c where k, b, and n are integer values\n");
}

void  FixedKBNApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "c:C:s:";

   AppendLongOpt(longOpts, "cmin",           required_argument, 0, 'c');
   AppendLongOpt(longOpts, "cmax",           required_argument, 0, 'C');
   AppendLongOpt(longOpts, "sequence",       required_argument, 0, 's');
}

parse_t FixedKBNApp::ParseOption(int opt, char *arg, const char *source)
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

      case 's':
         is_Sequence = arg;
         status = P_SUCCESS;
         break;
   }

   return status;
}

void FixedKBNApp::ValidateOptions(void)
{
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);

      iv_Terms.resize(il_MaxC - il_MinC + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);

      ProcessInputTermsFile(true);
   }
   else
   {
      if (il_MinC == 0)
         FatalError("cmin must be specified");

      if (il_MaxC == 0)
         FatalError("cmax must be specified");

      if (il_MaxC <= il_MinC)
         FatalError("cmax must be greater than cmin");

      if (is_Sequence.length() == 0)
         FatalError("sequence must be specified");

      if (sscanf(is_Sequence.c_str(), "%" SCNu64"*%u^%u+c", &il_K, &ii_Base, &ii_N) != 3)
         FatalError("sequence must be in form k*b^n+c where you specify values for k, b and n");

      if (ii_Base < 2)
         FatalError("base must be greater than 1");

      if (ii_N < 1)
         FatalError("n must be greater than 0");

      il_TermCount = il_MaxC - il_MinC + 1;

      iv_Terms.resize(il_MaxC - il_MinC + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);

      if ((il_K & 1) && (ii_Base & 1))
      {
         uint64_t c = il_MinC;
         
         // Remove all c that are odd if both k and b are odd since k*b^n+c will be divisible by 2.
         if (!(c & 1))
            c++;
         
         for (; c<=il_MaxC; c+=2)
         {
            iv_Terms[BIT(c)] = false;
            il_TermCount--;
         }
      }
      else
      {
         uint64_t c = il_MinC;
         
         // Remove all c that are even if both k and b are even since k*b^n+c will be divisible by 2.
         if (c & 1)
            c++;
         
         for (; c<=il_MaxC; c+=2)
         {
            iv_Terms[BIT(c)] = false;
            il_TermCount--;
         }
      }
   }

   char  fileName[50];

   if (is_OutputTermsFileName.length() == 0)
   {
      snprintf(fileName, sizeof(fileName), "k%" PRIu64"_b%u_n%u.pfgw", il_K,ii_Base, ii_N);
      is_OutputTermsFileName = fileName;
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

Worker *FixedKBNApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that FixedKBN inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new FixedKBNWorker(id, this);

   return theWorker;
}

void FixedKBNApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t bit;
   int64_t  c;
   uint64_t diff, lastPrime;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());

   if (!haveBitMap)
      il_MinC = il_MaxC = 0;

   if (!memcmp(buffer, "ABCD ", 5))
   {
      if (sscanf(buffer, "ABCD %" SCNu64"*%u^%d+$a  [%" SCNd64"] // Sieved to %" SCNu64"", &il_K, &ii_Base, &ii_N, &c, &lastPrime) != 5)
         FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());

      SetMinPrime(lastPrime);

      if (haveBitMap)
      {
         iv_Terms[c-il_MinC] = true;
         il_TermCount++;
      }
      else
         il_MinC = il_MaxC = c;
   }
   else
      FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (sscanf(buffer, "%" SCNu64 , &diff) != 1)
         FatalError("Line %s is malformed", buffer);

      c += diff;

      if (haveBitMap)
      {
         bit = BIT(c);
         iv_Terms[bit] = true;
         il_TermCount++;
      }
      else
      {
         if (il_MinC > c) il_MinC = c;
         if (il_MaxC < c) il_MaxC = c;
      }
   }

   fclose(fPtr);
}

bool FixedKBNApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint64_t k;
   uint32_t b, n;
   int64_t  c;

   if (sscanf(term, "%" SCNu64"*%u^%u%" SCNd64"", &k, &b, &n, &c) != 4)
      FatalError("Could not parse term %s", term);

   if (b != ii_Base)
      FatalError("Expected base %u in factor but found base %u", ii_Base, b);

   if (n != ii_N)
      FatalError("Expected n %u in factor but found %d", ii_N, n);

   if (k != il_K)
      FatalError("Expected k %" PRIu64" in factor but found %" PRIu64"", il_K, k);

   if (c < il_MinC || c > il_MaxC)
      return false;

   VerifyFactor(theFactor, c);

   uint64_t bit = c - il_MinC;

   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_TermCount--;

      return true;
   }

   return false;
}

void FixedKBNApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;

   // Don't save until we reach this point.  The term count will not match because
   // we are not locking the counters for small primes.
   if (largestPrime < GetMaxPrimeForSingleWorker())
      return;

   ip_FactorAppLock->Lock();

   int64_t c, previousC;
   uint64_t bit;

   c = il_MinC;

   bit = BIT(c);
   for (; c<=il_MaxC; c++)
   {
      if (iv_Terms[bit])
         break;

      bit++;
   }

   if (c > il_MaxC) {
      WriteToConsole(COT_OTHER, "No terms remaining therefore no output file was generated");
      return;
   }

   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   fprintf(termsFile, "ABCD %" PRIu64"*%u^%u+$a [%" PRId64"] // Sieved to %" PRIu64"\n", il_K, ii_Base, ii_N, c, largestPrime);

   previousC = c;
   termsCounted = 1;
   c++;

   bit = BIT(c);
   for (; c<=il_MaxC; c++)
   {
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%+" PRId64"\n", c - previousC);
         previousC = c;
         termsCounted++;
      }

      bit++;
   }

   fclose(termsFile);

   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void  FixedKBNApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   snprintf(extraText, maxTextLength, "%" PRId64" <= c <= %" PRId64", %" PRIu64 "*%u^%u+c", il_MinC, il_MaxC, il_K, ii_Base, ii_N);
}

bool  FixedKBNApp::ReportFactor(uint64_t theFactor, int64_t c, bool verifyFactor)
{
   bool     removedTerm = false;

   if (c < il_MinC || c > il_MaxC)
      return false;

   if (verifyFactor)
      VerifyFactor(theFactor, c);

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   int64_t bit = BIT(c);

   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      removedTerm = true;

      LogFactor(theFactor, "%" PRIu64"*%u^%u%+" PRId64"", il_K, ii_Base, ii_N, c);

      il_FactorCount++;
      il_TermCount--;
   }

   if (theFactor > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();

   return removedTerm;
}

void  FixedKBNApp::VerifyFactor(uint64_t theFactor, int64_t c)
{
   MpArith mp(theFactor);

   MpRes rem = mp.pow(mp.nToRes(ii_Base), ii_N);
   rem = mp.mul(rem, mp.nToRes(il_K));

   if (c > 0)
      rem = mp.add(rem, mp.nToRes((uint64_t)c));
   else if (c < 0)
      rem = mp.sub(rem, mp.nToRes((uint64_t)(-c)));

   if (rem != mp.zero())
      FatalError("%" PRIu64"*%u^%u%+" PRId64" mod %" PRIu64" = %" PRIu64"", il_K, ii_Base, ii_N, c, theFactor, mp.resToN(rem));
}
