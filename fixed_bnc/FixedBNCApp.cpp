/* FixedBNCApp.cpp -- (C) Mark Rodenkirch, November 2017

   Sieve for k*b^n+c for a range of k and fixed b, n, and c

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "FixedBNCApp.h"
#include "FixedBNCWorker.h"

#define APP_NAME        "fbncsieve"
#define APP_VERSION     "1.4"

#define BIT(k)          ((k) - il_MinK)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new FixedBNCApp();
}

FixedBNCApp::FixedBNCApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of k*b^n+c numbers for fixed b, n, and c and variable k");
   SetLogFileName("fbncsieve.log");

   is_Sequence = "";
   
   il_MinK = 0;
   il_MaxK = 0;
   it_Format = FF_ABCD;
   ib_Remove = false;
   
   il_MaxPrimeForValidFactor = PMAX_MAX_62BIT;
}

void FixedBNCApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-s --sequence=s       Sequence to find factors of in form k*b^n+c where b, n, and c are decimal values\n");
   printf("-f --format=f         Format of output file (A=ABC, D=ABCD (default), N=NEWPGEN)\n");
   printf("-r --remove           Remove k where k %% base = 0\n");
}

void  FixedBNCApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "rk:K:s:f:";

   AppendLongOpt(longOpts, "kmin",           required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",           required_argument, 0, 'K');
   AppendLongOpt(longOpts, "sequence",       required_argument, 0, 's');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
   AppendLongOpt(longOpts, "remove",         no_argument, 0, 'r');
}

parse_t FixedBNCApp::ParseOption(int opt, char *arg, const char *source)
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

      case 'f':
         char value;
         status = Parser::Parse(arg, "ADN", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'A')
            it_Format = FF_ABC;
         if (value == 'D')
            it_Format = FF_ABCD;
         if (value == 'N')
            it_Format = FF_NEWPGEN;
         break;
         
      case 's':
         is_Sequence = arg;
         status = P_SUCCESS;
         break;
         
      case 'r':
         ib_Remove = true;
         status = P_SUCCESS;
         break;
   }

   return status;
}

void FixedBNCApp::ValidateOptions(void)
{
   if (it_Format == FF_UNKNOWN)
      FatalError("File format not valid, use A (ABC), D (ABCD) or N (NewPGen)");
   
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      iv_Terms.resize(il_MaxK - il_MinK + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);
      
      ProcessInputTermsFile(true);

      if (ib_Remove)
      {
         uint64_t k = il_MinK;
         
         while (k % ii_Base > 0)
            k++;
         
         // Remvoe k that are divisible by the base
         for ( ; k<=il_MaxK; k+=ii_Base)
         {
            if (iv_Terms[BIT(k)])
            {
               il_TermCount--;
               iv_Terms[BIT(k)] = false;
            }
         }
      }
   }
   else
   {
      if (il_MinK == 0)
         FatalError("kmin must be specified");

      if (il_MaxK == 0)
         FatalError("kmax must be specified");
      
      if (il_MaxK <= il_MinK)
         FatalError("kmax must be greater than kmin");
      
      if (is_Sequence.length() == 0)
         FatalError("sequence must be specified");

      if (sscanf(is_Sequence.c_str(), "k*%u^%u%d", &ii_Base, &ii_N, &ii_C) != 3)
         FatalError("sequence must be in form k*b^n+c where you specify values for b, n and c");

      if (ii_Base < 2)
         FatalError("base must be greater than 1");

      if (ii_N < 1)
         FatalError("n must be greater than 0");
      
      if (ii_C != 1 && ii_C != -1)
         FatalError("c must be -1 or +1");
      
      il_TermCount = il_MaxK - il_MinK + 1;
      
      iv_Terms.resize(il_MaxK - il_MinK + 1);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);
      
      if (ib_Remove)
      {
         uint64_t k = il_MinK;
         
         while (k % ii_Base > 0)
            k++;
         
         // Remvoe k that are divisible by the base
         for ( ; k<=il_MaxK; k+=ii_Base)
         {
            il_TermCount--;
            iv_Terms[BIT(k)] = false;
         }
      }
   }
   
   if (is_OutputTermsFileName.length() == 0)
   {
      char  fileName[30];
      
      if (it_Format == FF_NEWPGEN)
         sprintf(fileName, "k_b%u_n%u%+d.npg", ii_Base, ii_N, ii_C);
      else
         sprintf(fileName, "k_b%u_n%u%+d.pfgw", ii_Base, ii_N, ii_C);
      
      is_OutputTermsFileName = fileName;
   }
   
   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   AdjustMaxPrime();
   
   // Allow only one worker to do work when processing the first chunk of
   // primes.  This will set the ib_BlockingForFirstChunk flag to true while
   // the first chunk of primes is being worked on and be set to false when 
   // that chunk is done.  This allows us to avoid locking when factors are 
   // reported, which significantly hurts performance as most terms will be
   // removed due to small primes.
   if (il_MinPrime < 1000)
      SetBlockWhenProcessingFirstChunk(true);
}

Worker *FixedBNCApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that FixedBNC inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new FixedBNCWorker(id, this);

   return theWorker;
}

void FixedBNCApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000], sign;
   uint32_t n, type;
   uint64_t bit, k, diff, lastPrime;
   format_t format = FF_UNKNOWN;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
   if (!haveBitMap)
      il_MinK = il_MaxK = 0;
   
   if (!memcmp(buffer, "ABCD ", 5))
   {
      if (sscanf(buffer, "ABCD $a*%u^%d%d  [%" SCNu64"] // Sieved to %" SCNu64"", &ii_Base, &ii_N, &ii_C, &k, &lastPrime) != 5)
         FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());
      
      format = FF_ABCD;

      if (haveBitMap)
      {
         iv_Terms[k-il_MinK] = true;
         il_TermCount++;
      }
      else
         il_MinK = il_MaxK = k;
   }
   else if (!memcmp(buffer, "ABC ", 4))
   {
      if (sscanf(buffer, "ABC $a*%u^%d%d // Sieved to %" SCNu64"", &ii_Base, &ii_N, &ii_C, &lastPrime) != 4)
         FatalError("Line 1 is not a valid ABC line in input file %s", is_InputTermsFileName.c_str());
      
      format = FF_ABC;
   }
   else if (buffer[0] > '0' && buffer[0] <= '9')
   {
      // handle newpgnen noise created by srsieve/srfile
      if ((sscanf(buffer, "%" SCNu64":%c:0:%u:%u", &lastPrime, &sign, &ii_Base, &type) != 4) &&
          (sscanf(buffer, "%" SCNu64":%c:1:%u:%u", &lastPrime, &sign, &ii_Base, &type) != 4))
         FatalError("Line 1 is not a valid newgpen line in input file %s", is_InputTermsFileName.c_str());
      
      ii_N = 0;
      ii_C = 0;
      
      if (sign == 'P' && type == 1)
         ii_C = +1;
      
      if (sign == 'M' && type == 2)
         ii_C = -1;
      
      if (ii_C == 0)
         FatalError("Line 1 is not a valid newgpen line in input file %s", is_InputTermsFileName.c_str());
      
      format = FF_NEWPGEN;
   }
   else
      FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
   
   SetMinPrime(lastPrime);
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;
   
      switch (format)
      {
         case FF_ABCD:
            if (sscanf(buffer, "%" SCNu64 , &diff) != 1)
               FatalError("Line %s is malformed", buffer);
            
            k += diff;
            break;

         case FF_ABC:
            if (sscanf(buffer, "%" SCNu64, &k) != 1)
               FatalError("Line %s is malformed", buffer);
            
            if (il_MinK == 0)
               il_MinK = il_MaxK = k;
            break;

         case FF_NEWPGEN:
            if (sscanf(buffer, "%" SCNu64 " %u", &k, &n) != 2)
               FatalError("Line %s is malformed", buffer);
            
            if (il_MinK == 0)
               il_MinK = il_MaxK = k;
            
            if (ii_N == 0)
               ii_N = n;
            
            if (ii_N != n)               
               FatalError("Input file %s has multiple values for k", is_InputTermsFileName.c_str());
            break;

         default:
            FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
      }
            
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

   fclose(fPtr);
}

bool FixedBNCApp::ApplyFactor(const char *term)
{
   uint64_t k;
   uint32_t b, n;
   int32_t  c;
   
   if (sscanf(term, "%" SCNu64"*%u^%u%u", &k, &b, &n, &c) != 4)
      FatalError("Could not parse term %s", term);

   if (b != ii_Base)
      FatalError("Expected base %u in factor but found base %u", ii_Base, b);
   
   if (n != ii_N)
      FatalError("Expected n %u in factor but found %d", ii_N, n);
   
   if (c != ii_C)
      FatalError("Expected c %+d in factor but found %d", ii_C, c);
        
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

void FixedBNCApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t kCount = 0;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
      
   if (it_Format == FF_ABCD)
      kCount = WriteABCDTermsFile(largestPrime, termsFile);
   
   if (it_Format == FF_ABC)
      kCount = WriteABCTermsFile(largestPrime, termsFile);
   
   if (it_Format == FF_NEWPGEN)
      kCount = WriteNewPGenTermsFile(largestPrime, termsFile);
   
   fclose(termsFile);
   
   if (kCount != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", kCount, il_TermCount);
   
   ip_FactorAppLock->Release();
}

uint64_t FixedBNCApp::WriteABCDTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0, previousK;
   uint64_t bit;

   k = il_MinK;
   
   bit = BIT(k);
   for (; k<=il_MaxK; k++)
   {      
      if (iv_Terms[bit])
         break;
      
      bit++;
   }

   if (k > il_MaxK)
      return 0;
   
   fprintf(termsFile, "ABCD $a*%u^%d%+d [%" PRIu64"] // Sieved to %" PRIu64"\n", ii_Base, ii_N, ii_C, k, maxPrime);
   
   previousK = k;
   kCount = 1;
   k++;
   
   bit = BIT(k);
   for (; k<=il_MaxK; k++)
   {
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%" PRIu64"\n", k - previousK);
         previousK = k;
         kCount++;
      }
      
      bit++;
   }

   return kCount;
}

uint64_t FixedBNCApp::WriteABCTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0;
   uint64_t bit;

   fprintf(termsFile, "ABC $a*%u^%u%+d // Sieved to %" PRIu64"\n", ii_Base, ii_N, ii_C, maxPrime);
      
   k = il_MinK;
   bit = BIT(k);

   for ( ; k<=il_MaxK; k++)
   {
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%" PRIu64"\n", k);
         kCount++;
      }
      
      bit++;
   }
   
   return kCount;
}

uint64_t FixedBNCApp::WriteNewPGenTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0;
   uint64_t bit;

   fprintf(termsFile, "%" PRIu64":%c:1:%u:%u\n", maxPrime, (ii_C == 1 ? 'P' : 'M'), ii_Base, (ii_C == 1 ? 1 : 2));
      
   k = il_MinK;
   bit = BIT(k);
   
   for ( ; k<=il_MaxK; k++)
   {
      if (iv_Terms[bit])
      {
         fprintf(termsFile, "%" PRIu64" %u\n", k, ii_N);
         kCount++;
      }
      
      bit++;
   }
   
   return kCount;
}

void  FixedBNCApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", k*%u^%u%+d", il_MinK, il_MaxK, ii_Base, ii_N, ii_C);
}

bool  FixedBNCApp::ReportFactor(uint64_t p, uint64_t k)
{
   bool     removedTerm = false;

   if (p > il_MaxPrimeForValidFactor)
      return false;
   
   if (!ib_BlockingForFirstChunk)
      ip_FactorAppLock->Lock();

   uint64_t bit = BIT(k);
      
   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      removedTerm = true;
      
      LogFactor(p, "%" PRIu64"*%u^%u%+d", k, ii_Base, ii_N, ii_C);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (!ib_BlockingForFirstChunk)
      ip_FactorAppLock->Release();
   
   return removedTerm;
}

// Don't sieve beyond sqrt(maxk*b^n+c).  The worker will not
// remove terms that are prime, so this means that all remaining
// terms will be prime when we reach p = sqrt(maxk*b^n+c).
void  FixedBNCApp::AdjustMaxPrime(void)
{  
   double    toobig = (double) GetMaxPrime();
   double    c = (double) ii_C;
   double    b = (double) ii_Base;
   double    bpown = 1.0;
   double    mink = (double) il_MinK;
   double    maxk = (double) il_MaxK;
   
   toobig *= toobig;
   
   for (uint32_t i=0; i<ii_N; i++)
   {
      bpown *= b;

      // If mink*b^n+1 > 2^128, then we don't need to worry
      // about sieving too deeply.
      if (mink*bpown+c > toobig)
         return;
   }

   double sqrtmax = sqrt(maxk*bpown+c);

   if (sqrtmax < toobig)
   {
      uint64_t maxp = (uint64_t) sqrtmax;
      
      // To addess truncating when converting from a double
      maxp++;
      
      // if maxp > sqrt(maxk*b^n+c), then stop sieving at sqrt(maxk*b^n+c)
      // since all remaining terms will be prime.
      if (maxp < GetMaxPrime())
      {
         if (maxp & 1)
            maxp += 1;
         
         SetMaxPrime(maxp, "All remaining terms will be prime");
         il_MaxPrimeForValidFactor = maxp;
      }
   }
}