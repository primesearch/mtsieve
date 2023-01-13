/* CunninghamChainApp.cpp -- (C) Mark Rodenkirch, June 2022

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
#include "../sieve/primesieve.hpp"
#include "CunninghamChainApp.h"
#include "CunninghamChainWorker.h"

#define APP_NAME        "ccsieve"
#define APP_VERSION     "1.1"

#define MAX_LENGTH      30
#define NMAX_MAX        (1 << 31)
#define MAX_FILES       9999

#define BIT_HK(k)       (((k) - il_MinK) >> 1)
#define BIT_AK(k)       ((k) - il_MinK)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new CunninghamChainApp();
}

CunninghamChainApp::CunninghamChainApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to eliminate terms for Cunningham Chain prime searches");
   SetLogFileName("ccsieve.log");
   
   it_Format = FF_CC;
   
   it_ChainKind = CCT_UNKNOWN;
   it_TermType = TT_UNKNOWN;
         
   ii_Base = 0;
   il_MinK = 0;
   il_MaxK = 0;
   ii_N = 0;
   ii_ChainLength = 0;
   il_Terms = 0;
   ii_Primes = 0;
   ib_HalfK = false;
   ii_NumberOfFiles = 1;
 
   SetAppMinPrime(3);
   
   iv_Terms.clear();
}

CunninghamChainApp::~CunninghamChainApp()
{
   if (ii_Primes != 0)
      xfree(ii_Primes);
   
   if (il_Terms != 0)
      xfree(il_Terms);
}

void CunninghamChainApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-c --chainkind=c      1 = first kind (first term end with -1), 2 = second kind (first term end with +1)\n");
   printf("-t --termtype=t       1 = b^n, 2 = primorial, 3 = factorial\n");
   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-l --length=l         Cunningham chain length\n");
   printf("-b --base=b           Base to search\n");
   printf("-n --n=n              n of b^n, n# for primorial, n! for factorial\n");
   printf("-f --format=f         Format of output file (C=CC (default), N=NEWPGEN)\n");
   printf("-N --numberOfFiles=N  Number of files to split k across\n");
}

void  CunninghamChainApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "b:k:K:m:n:f:c:t:N:l:";

   AppendLongOpt(longOpts, "chainkind",      required_argument, 0, 'c');
   AppendLongOpt(longOpts, "kmin",           required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",           required_argument, 0, 'K');
   AppendLongOpt(longOpts, "length",         required_argument, 0, 'l');
   AppendLongOpt(longOpts, "termtype",       required_argument, 0, 't');
   AppendLongOpt(longOpts, "base",           required_argument, 0, 'b');
   AppendLongOpt(longOpts, "n",              required_argument, 0, 'n');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
   AppendLongOpt(longOpts, "numberOfFiles",  required_argument, 0, 'N');
}

parse_t CunninghamChainApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;
   char value;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'c':
         status = Parser::Parse(arg, "12", value);
         
         it_ChainKind = CCT_UNKNOWN;
   
         if (value == '1')
            it_ChainKind = CCT_FIRSTKIND;
         if (value == '2')
            it_ChainKind = CCT_SECONDKIND;
         break;
      
      case 't':
         status = Parser::Parse(arg, "123", value);
         
         it_TermType = TT_UNKNOWN;
   
         if (value == '1')
            it_TermType = TT_BN;
         if (value == '2')
            it_TermType = TT_PRIMORIAL;
         if (value == '3')
            it_TermType = TT_FACTORIAL;
         break;
         
      case 'k':
         status = Parser::Parse(arg, 1, KMAX_MAX, il_MinK);
         break;

      case 'K':
         status = Parser::Parse(arg, 1, KMAX_MAX, il_MaxK);
         break;
         
      case 'l':
         status = Parser::Parse(arg, 2, MAX_LENGTH, ii_ChainLength);
         break;
         
      case 'b':
         status = Parser::Parse(arg, 2, 1 << 15, ii_Base);
         break;
         
      case 'n':
         status = Parser::Parse(arg, 3, NMAX_MAX, ii_N);
         break;

      case 'N':
         status = Parser::Parse(arg, 1, MAX_FILES, ii_NumberOfFiles);
         break;

      case 'f':
         status = Parser::Parse(arg, "CN", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'C')
            it_Format = FF_CC;
         if (value == 'N')
            it_Format = FF_NEWPGEN;
         break;
   }

   return status;
}

void CunninghamChainApp::ValidateOptions(void)
{
   if (it_Format == FF_UNKNOWN)
      FatalError("File format not valid, use A (ABC) or N (NewPGen)");
   
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);

      if (it_TermType == TT_PRIMORIAL)
         BuildPrimorialTerms();
      
      if (it_TermType == TT_FACTORIAL)
         BuildFactorialTerms();

      if (ib_HalfK)
         iv_Terms.resize(((il_MaxK - il_MinK) >> 1) + 1);
      else
         iv_Terms.resize((il_MaxK - il_MinK) + 1);
      
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);

      ProcessInputTermsFile(true);
   }
   else
   {
      if (it_ChainKind == CCT_UNKNOWN)
         FatalError("Chain kind must be specified");
      
      if (it_TermType == TT_UNKNOWN)
         FatalError("Term type must be specified");
      
      if (ii_ChainLength == 0)
         FatalError("Length must be specified");
      
      if (il_MinK == 0)
         FatalError("kmin must be specified");

      if (il_MaxK == 0)
         FatalError("kmax must be specified");
      
      if (il_MaxK <= il_MinK)
         FatalError("kmax must be greater than kmin");

      if (it_TermType == TT_BN && ii_Base == 0)
         FatalError("base must be specified");

      if (it_Format == FF_NEWPGEN && it_TermType == TT_FACTORIAL)
         FatalError("newpgen format not supported for factorials");
      
      if (il_MaxK < ii_N)
         FatalError("n must be less then maxk");
      
      il_TermCount = (il_MaxK - il_MinK) + 1;
         
      if (it_TermType == TT_BN)
      {
         if (il_MaxK < ii_Base)
            FatalError("base must be less then maxk");
         
         if (ii_Base == 2)
         {
            // We only care about odd k
            ib_HalfK = true;
            il_TermCount = (il_MaxK - il_MinK)/2 + 1;
         
            // Make minK odd
            if (!(il_MinK & 1))
               il_MinK++;
            
            // Make maxK odd
            if (!(il_MaxK & 1))
               il_MaxK--;
         }
         
         if (ii_Base & 0x01)
         {
            // We only care about even k
            ib_HalfK = true;
            il_TermCount = (il_MaxK - il_MinK)/2 + 1;
         
            // Make minK even
            if (il_MinK & 1)
               il_MinK++;
            
            // Make maxK even
            if (il_MaxK & 1)
               il_MaxK--;
         }
      }
      else
      {            
         if (it_TermType == TT_PRIMORIAL)
         {
            BuildPrimorialTerms();
            SetAppMinPrime(ii_N+1);
         }
         
         if (it_TermType == TT_FACTORIAL)
         {
            BuildFactorialTerms();
            SetAppMinPrime(ii_N+1);
         }
      }
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), true);      
   }

   if (is_OutputTermsFileName.length() == 0)
   {
      char filePrefix[60];
      
      if (it_TermType == TT_BN)
         sprintf(filePrefix, "cc_%u_%u", ii_Base, ii_N);
      
      if (it_TermType == TT_PRIMORIAL)
         sprintf(filePrefix, "cc_%up", ii_N);
      
      if (it_TermType == TT_FACTORIAL)
         sprintf(filePrefix, "cc_%uf", ii_N);

      is_OutputTermsFileName = filePrefix;
   }

   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(100000);
}

Worker *CunninghamChainApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that CunninghamChain inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new CunninghamChainWorker(id, this);

   return theWorker;
}

void CunninghamChainApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr;
   char     fileNamePattern[200], fileName[200];
   uint32_t fileCount = 0;

   il_TermCount = 0;

   if (!haveBitMap)
   {
      ii_Base = 0;
      ii_N = 0;
      ii_ChainLength = 0;
      il_MinK = PMAX_MAX_62BIT;
      il_MaxK = 0;
   }

   sprintf(fileName, "%s", is_InputTermsFileName.c_str());
      
   fileNamePattern[0] = 0;
   
   // Look for the pattern used when we have multiple output files
   for (uint32_t idx=0; idx<strlen(fileName); idx++)
   {
      if (fileName[idx] == '.' && fileName[idx+5] == '.')
      {
         fileName[idx] = 0;
         fileName[idx+5] = 0;
         sprintf(fileNamePattern, "%s.%s.%s", fileName, "%04u", &fileName[idx+6]);
         break;
      }
   }

   // If the file name has that pattern, then read from that set of files
   if (fileNamePattern[0] != 0)
   {
      for (uint32_t idx=1; idx<=MAX_FILES; idx++)
      {
         sprintf(fileName, fileNamePattern, idx);
         
         fPtr = fopen(fileName, "r");
      
         if (!fPtr)
            break;
         
         fileCount++;
         
         ProcessInputTermsFile(haveBitMap, fPtr, fileName, (fileCount == 1));
         
         fclose(fPtr);
      }
   }
   else
   {
      fPtr = fopen(fileName, "r");

      if (fPtr)
      {
         ProcessInputTermsFile(haveBitMap, fPtr, fileName, true);
         
         fileCount = 1;
         
         fclose(fPtr);
         
         if (haveBitMap)
            WriteToConsole(COT_OTHER, "Read input terms from 1 file");
         
         return;
      }
   }
   
   if (fileCount == 0)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (haveBitMap)
      WriteToConsole(COT_OTHER, "Read input terms from %d files", fileCount);
}

void CunninghamChainApp::ProcessInputTermsFile(bool haveBitMap, FILE *fPtr, char *fileName, bool firstFile)
{
   char       buffer[1000], *pos;
   int32_t    c = 2;
   uint32_t   base = 0, n = 0, chainLength;
   chainkind_t chainKind;
   termtype_t  termType = TT_UNKNOWN;
   uint64_t   k, lastPrime;
   format_t   format = FF_UNKNOWN;

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());

   pos = strstr(buffer, " //");
   if (pos)
   {
      *pos = 0;
      pos++;
      if (sscanf(pos, "// Sieved to %" SCNu64"", &lastPrime) == 1)
         SetMinPrime(lastPrime);
   }
   
   if (sscanf(buffer, "CC %u,%u,$a*%u^%u%d", &chainKind, &chainLength, &base, &n, &c) == 5)
   {
      format = FF_CC;
      termType = TT_BN;
   }
   else if (sscanf(buffer, "CC %u,%u,$a*%u#%d", &chainKind, &chainLength, &n, &c) == 4)
   {
      format = FF_CC;
      termType = TT_PRIMORIAL;
   }
   else if (sscanf(buffer, "CC %u,%u,$a*%u!%d", &chainKind, &chainLength, &ii_N, &c) == 4)
   {
      format = FF_CC;
      termType = TT_FACTORIAL;
   }
   else if (sscanf(buffer, "%" SCNu64":S:0:%u:74", &lastPrime, &base) == 2)
   {
      // The n is on each line
      format = FF_NEWPGEN;
      SetMinPrime(lastPrime);
      termType = TT_BN;
      chainLength = 2;
      chainKind = CCT_FIRSTKIND;
   }
   else if (sscanf(buffer, "%" SCNu64":C:0:%u:69", &lastPrime, &base) == 2)
   {
      // The n is on each line
      format = FF_NEWPGEN;
      SetMinPrime(lastPrime);
      termType = TT_BN;
      chainLength = 2;
      chainKind = CCT_SECONDKIND;
   }
   else if (sscanf(buffer, "%" SCNu64":1:%u:%u:1066", &lastPrime, &chainLength, &base) == 3)
   {
      format = FF_NEWPGEN;
      SetMinPrime(lastPrime);
      termType = TT_BN;
      chainKind = CCT_FIRSTKIND;
   }
   else if (sscanf(buffer, "%" SCNu64":1:%u:%u:1045", &lastPrime, &chainLength, &base) == 3)
   {
      format = FF_NEWPGEN;
      SetMinPrime(lastPrime);
      termType = TT_BN;
      chainKind = CCT_SECONDKIND;
   }
   else if (sscanf(buffer, "%" SCNu64":1:%u:2:106", &lastPrime, &chainLength) == 2)
   {
      // The primorial is on each line
      format = FF_NEWPGEN;
      termType = TT_PRIMORIAL;
      
      if (firstFile)
      {
         SetMinPrime(lastPrime);
         chainKind = CCT_FIRSTKIND;
      }
   }
   else if (sscanf(buffer, "%" SCNu64":2:%u:2:85", &lastPrime, &chainLength) == 2)
   {
      // The primorial is on each line
      format = FF_NEWPGEN;
      termType = TT_PRIMORIAL;
      
      if (firstFile)
      {
         SetMinPrime(lastPrime);
         chainKind = CCT_SECONDKIND;
      }
   }
   else
      FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());    
   
   if (firstFile)
   {
      it_ChainKind = chainKind;
      ii_ChainLength = chainLength;
      it_TermType = termType;
      ii_Base = base;
      ii_N = n;
   }
   
   if (format == FF_CC && chainKind == CCT_FIRSTKIND && c != -1)
      FatalError("Chains of the first kind must end with -1");
   
   if (format == FF_CC && chainKind == CCT_SECONDKIND && c != +1)
      FatalError("Chains of the first kind must end with +1");
   
   if (base != ii_Base)
      FatalError("Mixed bases in input files");
   
   if (format == FF_CC && n != ii_N)
      FatalError("Mixed n in input files");
      
   if (chainLength != ii_ChainLength)
      FatalError("Mixed length in input files");
   
   if (chainKind != it_ChainKind)
      FatalError("Mixed tyeps of chains in input files");
   
   if (termType != it_TermType)
      FatalError("Mixed types of terms in input files");

   // We only care about odd k or even k, but not both
   if (it_TermType == TT_BN && (ii_Base == 2 || ii_Base & 1))
      ib_HalfK = true;
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (format == FF_CC)
      {
         if (sscanf(buffer, "%" SCNu64"" , &k) != 1)
            FatalError("Line %s is malformed", buffer);
      }
      else
      {
         if (sscanf(buffer, "%" SCNu64" %u" , &k, &n) != 2)
            FatalError("Line %s is malformed", buffer);

         if (ii_N == 0)
            ii_N = n;
         
         if (ii_N != n)
            FatalError("Mixed n values in input file");
      }
            
      if (haveBitMap)
      {
         if (ib_HalfK)
            iv_Terms[BIT_HK(k)] = true;
         else
            iv_Terms[BIT_AK(k)] = true;
         
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

void CunninghamChainApp::BuildPrimorialTerms(void)
{
   primesieve::iterator   primeIterator;
   uint32_t termIdx = 0;
   
   primeIterator.skipto(1, ii_N);
   
   uint32_t primeCount = primesieve::count_primes(1, ii_N);
   ii_Primes = (uint32_t *) xmalloc((10 + primeCount) * sizeof(uint32_t));
   primeCount = 0;
   
   while (true)
   {
      ii_Primes[primeCount] = (uint32_t) primeIterator.next_prime();
      
      if (ii_Primes[primeCount] == ii_N)
         break;
      
      if (ii_Primes[primeCount] > ii_N)
         FatalError("%u is not prime, consider %u or %u\n", ii_N, ii_Primes[primeCount-1], ii_Primes[primeCount]);
         
      primeCount++;
   }
   
   il_Terms = (uint64_t *) xmalloc((10 + (primeCount / 3)) * sizeof(uint64_t));
   
   il_Terms[0] = 2;
   termIdx = 0;
   for (uint32_t idx=1; ii_Primes[idx] > 0; idx++)
   {
      uint64_t mult = il_Terms[termIdx] * ii_Primes[idx];
      
      // If we overflows then mult will be less than il_Terms[termIdx]
      if ((double) il_Terms[termIdx] * (double) ii_Primes[idx] < (double) PMAX_MAX_62BIT)
      {
         il_Terms[termIdx] = mult;
         continue;
      }
      
      termIdx++;
      il_Terms[termIdx] = ii_Primes[idx];
   }
}

void CunninghamChainApp::BuildFactorialTerms(void)
{
   primesieve::iterator   primeIterator;
   uint32_t termIdx = 0;

   il_Terms = (uint64_t *) xmalloc((10 + (ii_N / 3)) * sizeof(uint64_t));

   for (uint32_t n=2; n<=ii_N; n++)
   {
      if (il_Terms[termIdx] == 0)
      {
         il_Terms[termIdx] = n;
         continue;
      }
      
      uint64_t mult = il_Terms[termIdx] * n;
      
      // If we overflows then mult will be less than il_Terms[termIdx]
      if ((double) il_Terms[termIdx] * (double) n < (double) PMAX_MAX_62BIT)
      {
         il_Terms[termIdx] = mult;
         continue;
      }
      
      termIdx++;
      il_Terms[termIdx] = n;
   }
}

bool CunninghamChainApp::ApplyFactor(uint64_t theFactor, const char *term)
{   
   FatalError("ApplyFactor is not implemented yet");
   
   return false;
}

void CunninghamChainApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint64_t k = il_MinK;
   
   ip_FactorAppLock->Lock();
   
   if (il_TermCount == 0)
      FatalError("No remaining terms");

   if (ii_NumberOfFiles > 1 && il_TermCount < 10 * ii_NumberOfFiles)
   {
      WriteToConsole(COT_OTHER, "Reducing to 1 file as we need at least 10 terms per file");
      ii_NumberOfFiles = 1;
   }
   
   for (uint32_t idx=1; idx<=ii_NumberOfFiles; idx++)
   {
      char fileName[200];
      
      if (ii_NumberOfFiles == 1)
         sprintf(fileName, "%s.%s", is_OutputTermsFileName.c_str(), (it_Format == FF_CC ? "cc" : "npg"));
      else
         sprintf(fileName, "%s.%04u.%s", is_OutputTermsFileName.c_str(), idx, (it_Format == FF_CC ? "cc" : "npg"));
      
      FILE    *termsFile = fopen(fileName, "w");

      if (!termsFile)
         FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
      
      if (it_Format == FF_CC)
         termsCounted += WriteCCTermsFile(largestPrime, termsFile, k);
      
      if (it_Format == FF_NEWPGEN)
         termsCounted += WriteNewPGenTermsFile(largestPrime, termsFile, k);
      
      fclose(termsFile);
      
      if (k == 0)
         break;
   }
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

uint64_t CunninghamChainApp::WriteCCTermsFile(uint64_t largestPrime, FILE *termsFile, uint64_t &nextK)
{
   uint64_t k, kCount = 0;
   uint64_t kToWrite = 1 + (il_TermCount / ii_NumberOfFiles);
   char     term[50];

   if (it_TermType == TT_BN)
      sprintf(term, "$a*%u^%u%+d", ii_Base, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
   
   if (it_TermType == TT_PRIMORIAL)
      sprintf(term, "$a*%u#%+d", ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
   
   if (it_TermType == TT_FACTORIAL)
      sprintf(term, "$a*%u!%+d", ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
      
   fprintf(termsFile, "CC %u,%u,%s // Sieved to %" SCNu64"\n", it_ChainKind, ii_ChainLength, term, largestPrime);

   if (ib_HalfK)
   {
      for (k=nextK; k<=il_MaxK; k+=2)
      {
         if (iv_Terms[BIT_HK(k)])
         {
            if (kCount == kToWrite)
            {
               nextK = k;
               return kCount;
            }
         
            fprintf(termsFile, "%" PRIu64"\n", k);
            
            kCount++;
         }
      }
   }
   else
   {
      for (k=nextK; k<=il_MaxK; k++)
      {
         if (iv_Terms[BIT_AK(k)])
         {
            if (kCount == kToWrite)
            {
               nextK = k;
               return kCount;
            }
         
            fprintf(termsFile, "%" PRIu64"\n", k);
            
            kCount++;
         }
      }
   }

   nextK = 0;
   return kCount;
}

uint64_t CunninghamChainApp::WriteNewPGenTermsFile(uint64_t largestPrime, FILE *termsFile, uint64_t &nextK)
{
   uint64_t k, kCount = 0;
   uint64_t kToWrite = 1 + (il_TermCount / ii_NumberOfFiles);

   if (ii_Base > 0)
   {
      if (it_ChainKind == 1)
      {
         if (ii_ChainLength == 2)
            fprintf(termsFile, "%" PRIu64":S:0:%u:74", largestPrime, ii_Base);
         else
            fprintf(termsFile, "%" PRIu64":1:%u:%u:1066", largestPrime, ii_ChainLength, ii_Base);
      }
      
      if (it_ChainKind == 2)
      {
         if (ii_ChainLength == 2)
            fprintf(termsFile, "%" PRIu64":C:0:%u:69", largestPrime, ii_Base);
         else
            fprintf(termsFile, "%" PRIu64":2:%u:%u:1045", largestPrime, ii_ChainLength, ii_Base);
      }
   }
   else
   {
      if (it_ChainKind == 1)
         fprintf(termsFile, "%" PRIu64":1:%u:2:106", largestPrime, ii_ChainLength);
      
      if (it_ChainKind == 2)
         fprintf(termsFile, "%" PRIu64":2:%u:2:85", largestPrime, ii_ChainLength);
   }   

   if (ib_HalfK)
   {
      for (k=il_MinK; k<=il_MaxK; k+=2)
      {
         if (iv_Terms[BIT_HK(k)])
         {
            if (kCount == kToWrite)
            {
               nextK = k;
               return kCount;
            }
         
            fprintf(termsFile, "%" PRIu64" %u\n", k, ii_N);
            
            kCount++;
         }
      }
   }
   else
   {
      for (k=il_MinK; k<=il_MaxK; k++)
      {
         if (iv_Terms[BIT_AK(k)])
         {
            if (kCount == kToWrite)
            {
               nextK = k;
               return kCount;
            }
         
            fprintf(termsFile, "%" PRIu64" %u\n", k, ii_N);
            
            kCount++;
         }
      }
   }
   
   nextK = 0;
   return kCount;
}

void  CunninghamChainApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   int c = (it_ChainKind == 1 ? -1 : +1);
   
   if (it_TermType == TT_BN)
      sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", first term of k*%u^%u%+d for CC length %u", il_MinK, il_MaxK, ii_Base, ii_N, c, ii_ChainLength);
   
   if (it_TermType == TT_PRIMORIAL)
      sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", first term of k*%u#%+d for CC length %u", il_MinK, il_MaxK, ii_N, c, ii_ChainLength);
   
   if (it_TermType == TT_FACTORIAL)
      sprintf(extraTtext, "%" PRIu64 " < k < %" PRIu64", first term of k*%u!%+d for CC length %u", il_MinK, il_MaxK, ii_N, c, ii_ChainLength);
}

void  CunninghamChainApp::ReportFactor(uint64_t theFactor, uint64_t k, uint32_t termInChain)
{
   // If the first term is valid, then the rest are valid.  In other words 
   // k*x (mod p) = (k+p)*x (mod p) = ... = (k+n*p)*x (mod p)
   VerifyFactor(theFactor, k, termInChain);
   
   ip_FactorAppLock->Lock();
   
   do
   {
      uint64_t bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));

      if (iv_Terms[bit])
      {
         iv_Terms[bit] = false;
         
         if (if_FactorFile != 0)
         {
            char     term[50];

            if (termInChain == 1)
            {
               if (it_TermType == TT_BN)
                  sprintf(term, "%" PRIu64"*%u^%u%+d", k, ii_Base, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
                  
               if (it_TermType == TT_PRIMORIAL)
                  sprintf(term, "%" PRIu64"*%u#%+d", k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
            
               if (it_TermType == TT_PRIMORIAL)
                  sprintf(term, "%" PRIu64"*%u#%+d", k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
            }
            else
            {
               uint32_t mult = 1 << (termInChain - 1);
               int32_t add = (it_ChainKind == CCT_FIRSTKIND ? (mult - 1) : -(mult - 1));
               
               if (it_TermType == TT_BN)
                  sprintf(term, "%u*(%" PRIu64"*%u^%u%+d)%+d", mult, k, ii_Base, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1), add);
                  
               if (it_TermType == TT_PRIMORIAL)
                  sprintf(term, "%u*(%" PRIu64"*%u#%+d)%+d", mult, k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1), add);
            
               if (it_TermType == TT_PRIMORIAL)
                  sprintf(term, "%u*(%" PRIu64"*%u#%+d)%+d", mult, k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1), add);
            }
            
            LogFactor(theFactor, "%s", term);
         }
         
         il_FactorCount++;
         il_TermCount--;
      }
      
      if (ib_HalfK)
      {
         // We only care about every other k
         k += (theFactor << 1);
      }
      else
         k += theFactor; 
   } while (k <= il_MaxK);
   
   ip_FactorAppLock->Release();
}

void  CunninghamChainApp::VerifyFactor(uint64_t theFactor, uint64_t k, uint32_t termInChain)
{
   MpArith  mp(theFactor);
   MpRes    pOne = mp.one();
   MpRes    resRem = pOne;
   
   if (it_TermType == TT_BN)
   {
      resRem = mp.pow(mp.nToRes(ii_Base), ii_N);
   }
   else
   {   
      uint32_t idx = 0;
      while (il_Terms[idx] > 0)
      {
         resRem = mp.mul(resRem, mp.nToRes(il_Terms[idx]));

         idx++;
      }
   }

   resRem = mp.mul(resRem, mp.nToRes(k));
   
   if (it_ChainKind == CCT_FIRSTKIND)
      resRem = mp.sub(resRem, pOne);
   else
      resRem = mp.add(resRem, pOne);

   if (termInChain > 1)
   {
      uint32_t mult = 1 << (termInChain - 1);
      
      resRem = mp.mul(resRem, mp.nToRes(mult));
      
      if (it_ChainKind == CCT_FIRSTKIND)
         resRem = mp.add(resRem, mp.nToRes(mult - 1));
      else
         resRem = mp.sub(resRem, mp.nToRes(mult - 1));
   }
         
   if (resRem != mp.zero())
   {
      char     term[50];

      if (termInChain == 1)
      {
         if (it_TermType == TT_BN)
            sprintf(term, "%" PRIu64"*%u^%u%+d", k, ii_Base, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
            
         if (it_TermType == TT_PRIMORIAL)
            sprintf(term, "%" PRIu64"*%u#%+d", k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
      
         if (it_TermType == TT_FACTORIAL)
            sprintf(term, "%" PRIu64"*%u!%+d", k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1));
      }
      else
      {
         uint32_t mult = 1 << (termInChain - 1);
         int32_t add = (it_ChainKind == CCT_FIRSTKIND ? (mult - 1) : -(mult - 1));
         
         if (it_TermType == TT_BN)
            sprintf(term, "%u*(%" PRIu64"*%u^%u%+d)%+d", mult, k, ii_Base, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1), add);
            
         if (it_TermType == TT_PRIMORIAL)
            sprintf(term, "%u*(%" PRIu64"*%u#%+d)%+d", mult, k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1), add);
      
         if (it_TermType == TT_FACTORIAL)
            sprintf(term, "%u*(%" PRIu64"*%u!%+d)%+d", mult, k, ii_N, (it_ChainKind == CCT_FIRSTKIND ? -1 : +1), add);
      }
      
      FatalError("Invalid factor: %" PRIu64" is not a factor of %s", theFactor, term);
   }
}
