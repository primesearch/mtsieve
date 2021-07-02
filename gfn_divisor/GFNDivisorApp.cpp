/* GFNDivisorApp.cpp -- (C) Mark Rodenkirch, November 2017

   Sieve for k*2^n+1 for a range of k and a range of n

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <time.h>
#include "../core/Parser.h"
#include "../core/MpArithVector.h"
#include "../sieve/primesieve.hpp"
#include "GFNDivisorApp.h"
#include "GFNDivisorWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

#define APP_VERSION     "2.0"

#ifdef HAVE_GPU_WORKERS
#include "GFNDivisorGpuWorker.h"
#define APP_NAME        "gfndsievecl"
#else
#define APP_NAME        "gfndsieve"
#endif

#define KMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

#define BIT(k)          (((k) - il_MinK) >> 1)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new GFNDivisorApp();
}

GFNDivisorApp::GFNDivisorApp(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of k*2^n+1 numbers for variable k and n");
   SetLogFileName("gfndsieve.log");

   il_MinKOriginal = il_MinK = 0;
   il_MaxKOriginal = il_MaxK = 0;
   ii_MinNOriginal = ii_MinN = 0;
   ii_MaxNOriginal = ii_MaxN = 0;
      
   ii_NsPerFile = NMAX_MAX;
   ib_TestTerms = false;
   il_KPerChunk = 10000000000L;
   ib_UseTermsBitmap = true;
   ii_SmallPrimeFactorLimit = 32767;

   SetAppMinPrime(3);
   
#ifdef HAVE_GPU_WORKERS
   ii_MaxGpuSteps = 1000000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}

void GFNDivisorApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k              minimum k to search\n");
   printf("-K --kmax=K              maximum k to search\n");
   printf("-n --nmin=N              minimum n to search\n");
   printf("-N --nmax=N              maximum n to search\n");
   printf("-T --nsperfile=T         number of n per output file\n");
   printf("-x --testterms           test remaining terms for GFN divisibility\n");
   printf("-X --termsperchunk=X     used with -x, number of terms to sieve at a time (default 1e10)\n");
   printf("-r --notermsbitmap       do not generate terms bitmap\n");
   printf("-R --smallprimelimit=R   used with -r, do not output terms with a divisor < R (default 32767)\n");
   
#ifdef HAVE_GPU_WORKERS
   printf("-S --step=S              max steps iterated per call to GPU (default %u)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M        max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  GFNDivisorApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:n:N:c:T:xX:rR:";

   AppendLongOpt(longOpts, "kmin",              required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",              required_argument, 0, 'K');
   AppendLongOpt(longOpts, "nmin",              required_argument, 0, 'n');
   AppendLongOpt(longOpts, "nmax",              required_argument, 0, 'N');
   AppendLongOpt(longOpts, "nsperfile",         required_argument, 0, 'T');
   AppendLongOpt(longOpts, "testterms",         no_argument,       0, 'x');
   AppendLongOpt(longOpts, "termsperchunk",     required_argument, 0, 'X');
   AppendLongOpt(longOpts, "notermsbitmap",     no_argument,       0, 'r');
   AppendLongOpt(longOpts, "smallprimelimit",   required_argument, 0, 'R');
   
#ifdef HAVE_GPU_WORKERS
   shortOpts += "S:M:";
   
   AppendLongOpt(longOpts, "maxsteps",       required_argument, 0, 'S');
   AppendLongOpt(longOpts, "maxfactors",     required_argument, 0, 'M');
#endif
}

parse_t GFNDivisorApp::ParseOption(int opt, char *arg, const char *source)
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
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MinN);
         break;

      case 'N':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_MaxN);
         break;

      case 'r':
         ib_UseTermsBitmap = false;
         status = P_SUCCESS;
         break;
         
      case 'R':
         status = Parser::Parse(arg, 1, 1000000, ii_SmallPrimeFactorLimit);
         break;
         
      case 'T':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_NsPerFile);
         break;
         
      case 'x':
         ib_TestTerms = true;
         status = P_SUCCESS;
         break;
         
      case 'X':
         status = Parser::Parse(arg, 1000000000, KMAX_MAX, il_KPerChunk);
         break;

#ifdef HAVE_GPU_WORKERS
      case 'S':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxGpuSteps);
         break;
         
      case 'M':
         status = Parser::Parse(arg, 1, 1000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void GFNDivisorApp::ValidateOptions(void)
{
   uint64_t kCount;
   uint32_t nCount;
   char     maxPrimeStr[50];
   
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "gfnd";

   // We need to ensure compatibility between GMP and the x86_asm_ext functions
   if (sizeof(uint64_t) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
  
#ifdef WIN32
   if (sizeof(unsigned long long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#else
   if (sizeof(unsigned long) != sizeof(mp_limb_t))
     FatalError("GMP limb size is not 64 bits");
#endif

   if (!ib_UseTermsBitmap)
   {
      WriteToConsole(COT_OTHER, "Disabling output of remaining terms due to no bitmap");
      
      iv_SmallPrimes.clear();
      
      // Both the CPU and GPU have special logic to avoid reporting
      // factors when the candidate has a factor < 50.
      uint32_t primesInRange = primesieve::count_primes(50, ii_SmallPrimeFactorLimit);
      
      // We want a number of primes that is divisble by 4.
      while (primesInRange & 0x03)
         primesInRange++;
      
      // Generate the list of small primes
      primesieve::generate_n_primes(primesInRange, 50, &iv_SmallPrimes);

      if (ib_TestTerms)
      {
         ib_TestTerms = false;
         WriteToConsole(COT_OTHER, "Disabling testing of terms due to no bitmap");
      }
      
      if (is_OutputFactorsFileName.length() == 0)
         is_OutputFactorsFileName = "gfnd_fact.txt";
   }

   is_OutputTermsFilePrefix = is_OutputTermsFileName;

   if (is_InputTermsFileName.length() > 0)
   {
      if (ib_TestTerms)
         FatalError("cannot use -i and -x together");
      
      ProcessInputTermsFile(false);
      
      // Make minK odd
      if (!(il_MinK & 1))
         il_MinK++;
      
      // Make maxK odd
      if (!(il_MaxK & 1))
         il_MaxK--;
      
      kCount = (il_MaxK - il_MinK)/2 + 1;
      nCount = ii_MaxN - ii_MinN + 1;
         
      if (ib_UseTermsBitmap)
      {
         iv_Terms.resize(nCount);
         
         for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
         {
            iv_Terms[n-ii_MinN].resize(kCount);
            std::fill(iv_Terms[n-ii_MinN].begin(), iv_Terms[n-ii_MinN].end(), false);
         }
               
         ProcessInputTermsFile(true);
      }
   }
   else
   {
      if (il_MinK == 0)
         FatalError("kmin must be specified");

      if (il_MaxK == 0)
         FatalError("kmax must be specified");
      
      if (ii_MinN == 0)
         FatalError("nmin must be specified");

      if (ii_MaxN == 0)
         FatalError("nmax must be specified");
      
      if (il_MaxK <= il_MinK)
         FatalError("kmax must be greater than kmin");
      
      if (ii_MaxN < ii_MinN)
         FatalError("nmax must be greater than nmin");
      
      // Make minK odd
      if (!(il_MinK & 1))
         il_MinK++;
      
      // Make maxK odd
      if (!(il_MaxK & 1))
         il_MaxK--;
      
      nCount = ii_MaxN - ii_MinN + 1;
   }
   
   if (nCount / ii_NsPerFile > 9999)
      FatalError("nsperfile is too small as too many files would be created");

   if (ib_TestTerms)   
      ip_GFNDivisorTester = new GFNDivisorTester(this);

   // Not expecting anyone to use gfndsieve for n > 2000 since gfndsieve+pfgw should be faster
   if (ib_TestTerms && il_MaxPrime == il_AppMaxPrime) {
      if (ii_MaxN <= 100000) il_MaxPrime = 30000000000;
      if (ii_MaxN <=  50000) il_MaxPrime = 30000000000;
      if (ii_MaxN <=  30000) il_MaxPrime = 30000000000;
      if (ii_MaxN <=  20000) il_MaxPrime = 30000000000;
      if (ii_MaxN <=  10000) il_MaxPrime = 30000000000;
      if (ii_MaxN <=   5000) il_MaxPrime = 30000000000;
      if (ii_MaxN <=   3000) il_MaxPrime = 10000000000;
      if (ii_MaxN <=   2000) il_MaxPrime = 10000000000;
      if (ii_MaxN <=   1000) il_MaxPrime = 10000000000;
      if (ii_MaxN <=    800) il_MaxPrime = 10000000000;
      if (ii_MaxN <=    600) il_MaxPrime = 10000000000;
      if (ii_MaxN <=    400) il_MaxPrime =  3000000000;
      if (ii_MaxN <=    300) il_MaxPrime =  3000000000;
      if (ii_MaxN <=    200) il_MaxPrime =  1000000000;
      if (ii_MaxN <=    150) il_MaxPrime =  1000000000;
      if (ii_MaxN <=    100) il_MaxPrime =  1000000000;
      if (ii_MaxN <=     80) il_MaxPrime =  1000000000;
      if (ii_MaxN <=     60) il_MaxPrime =   100000000;
      if (ii_MaxN <=     40) il_MaxPrime =     3000000;

      if (il_MaxPrime != il_AppMaxPrime)
      {
         SetMaxPrime(il_MaxPrime, "testing terms for GFN divisibility");
         
         sprintf(maxPrimeStr, "%" PRIu64"", il_MaxPrime);
         
         WriteToConsole(COT_OTHER, "Sieving to %s due to testing of terms", maxPrimeStr);
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
   
   // We want each p to divide a maximum of one k per n because the GPU kernel won't support more.
   SetMinGpuPrime((il_MaxK-il_MinK) >> 1);
}

void  GFNDivisorApp::PreSieveHook(void)
{
   uint64_t kCount;
   uint32_t nCount;
   
   kCount = (il_MaxK - il_MinK)/2 + 1;
   nCount = ii_MaxN - ii_MinN + 1;
   il_TotalTerms = kCount * nCount;
   
   if (!ib_UseTermsBitmap)
      return;
   
   if (is_InputTermsFileName.length() > 0)
      return;
   
   if (!ib_TestTerms)
   {      
      il_TermCount = il_TotalTerms;
      
      iv_Terms.resize(nCount);
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         iv_Terms[n-ii_MinN].resize(kCount);
         std::fill(iv_Terms[n-ii_MinN].begin(), iv_Terms[n-ii_MinN].end(), true);
      }

      return;
   }
         
   ii_NPerChunk = 1;
   
   if (il_KPerChunk > (il_MaxK - il_MinK)/2 + 1)
      ii_NPerChunk = il_KPerChunk / ((il_MaxK - il_MinK)/2 + 1);

   if (ii_NPerChunk > nCount)
      ii_NPerChunk = nCount;
   
   if (il_MinKOriginal == 0)
   {
      iv_Terms.resize(ii_NPerChunk);
      
      for (uint32_t n=0; n<ii_NPerChunk; n++)
         iv_Terms[n].resize(il_KPerChunk);
      
      il_MinKOriginal = il_MinK;
      il_MaxKOriginal = il_MaxK;
      il_MinKInChunk = il_MinK;
      
      ii_MinNOriginal = ii_MinN;
      ii_MaxNOriginal = ii_MaxN;
      ii_MinNInChunk = ii_MinN;
      
   }
   
   for (uint32_t n=0; n<ii_NPerChunk; n++)
      std::fill(iv_Terms[0].begin(), iv_Terms[0].end(), true);

   // We want il_MinK and il_MaxK to be set to the correct range
   // of k before we start sieving.
   il_MinK = il_MinKInChunk;
   il_MaxK = il_MinK + il_KPerChunk + 1;
   
   ii_MinN = ii_MinNInChunk;
   ii_MaxN = ii_MinN + ii_NPerChunk - 1;

   uint32_t nInChunk = ii_NPerChunk;
   uint64_t kInChunk = il_KPerChunk;
   
   if (il_MaxK > il_MaxKOriginal)
   {
      il_MaxK = il_MaxKOriginal;
      kInChunk = (il_MaxK - il_MinK)/2 + 1;
   }
   
   if (ii_MaxN > ii_MaxNOriginal)
   {
      ii_MaxN = ii_MaxNOriginal;
      nInChunk = ii_MaxN - ii_MinN + 1;
   }
   
   il_TotalTermsInChunk = il_TermCount = kInChunk * nInChunk;
   il_FactorCount = 0;
   
   ip_GFNDivisorTester->StartedSieving();
}

bool  GFNDivisorApp::PostSieveHook(void)
{
   if (!ib_TestTerms)
      return true;
   
   ip_GFNDivisorTester->TestRemainingTerms(il_TotalTerms, il_TotalTermsInChunk, il_TermCount);
   
   // Set the starting k for the next range to be sieved.
   il_MinKInChunk = il_MaxK;

   if (il_MinKInChunk < il_MaxKOriginal)
      return false;
   
   // If we finished all k for this n, then we are done if we finished
   // all n as well
   if (ii_MaxN >= ii_MaxNOriginal)
      return true;
   
   // Set k and n for the next range
   
   // I know that we will resieve a single k for each successive chunk.
   // This is done so the output looks nicer, i.e. 1e9 to 2e9 rather than
   // 1e9 to 1999999999 then 2e9 to 2999999999.
   il_MinKInChunk = il_MinKOriginal;
   ii_MinNInChunk = ii_MaxN + 1;
   
   return false;
}

Worker *GFNDivisorApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
#ifdef HAVE_GPU_WORKERS
   if (gpuWorker)
      return new GFNDivisorGpuWorker(id, this);
#endif

   return new GFNDivisorWorker(id, this);
}

void GFNDivisorApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr;
   char     fileName[200];
   uint32_t fileCount = 0;

   il_TermCount = 0;
      
   if (!haveBitMap)
   {
      ii_MinN = ii_MaxN = 0;
      il_MinK = il_MaxK = 0;
   }
   
   sprintf(fileName, "%s", is_InputTermsFileName.c_str());
   
   fPtr = fopen(fileName, "r");
   
   if (!fPtr)
   {
      sprintf(fileName, "%s.pfgw", is_InputTermsFileName.c_str());
      
      fPtr = fopen(fileName, "r");
   }

   if (fPtr)
   {
      ProcessInputTermsFile(haveBitMap, fPtr, fileName, true);
      
      fclose(fPtr);
      
      if (haveBitMap)
         WriteToConsole(COT_OTHER, "Read input terms from 1 file");
      
      return;
   }
      
   for (int i=1; i<10000; i++)
   {
      sprintf(fileName, "%s_%04d.pfgw", is_InputTermsFileName.c_str(), i);
      
      fPtr = fopen(fileName, "r");
   
      if (!fPtr)
         break;
      
      fileCount++;
      
      ProcessInputTermsFile(haveBitMap, fPtr, fileName, (fileCount == 1));
      
      fclose(fPtr);
   }
   
   if (fileCount == 0)
      FatalError("Unable to open any input file with the prefix %s", is_InputTermsFileName.c_str());

   if (haveBitMap)
      WriteToConsole(COT_OTHER, "Read input terms from %d files", fileCount);
}

void GFNDivisorApp::ProcessInputTermsFile(bool haveBitMap, FILE *fPtr, char *fileName, bool firstFile)
{
   char     buffer[1000];
   uint32_t n;
   uint64_t k, diff, minPrime;

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", fileName);
   
   if (!memcmp(buffer, "ABCD ", 5))
   {
      if (sscanf(buffer, "ABCD $a*2^%u+1 [%" SCNu64"] // Sieved to %" SCNu64"", &n, &k, &minPrime) != 3)
         FatalError("Line 1 is not a valid ABCD line in input file %s", fileName);
      
      SetMinPrime(minPrime);

      if (haveBitMap)
      {
         iv_Terms[n-ii_MinN][BIT(k)] = true;
         il_TermCount++;
      }
      else
      {
         if (firstFile)
         {
            il_MinK = il_MaxK = k;
            ii_MinN = ii_MaxN = n;
         }
         else
         {
            if (il_MinK > k) il_MinK = k;
            if (il_MaxK < k) il_MaxK = k;
            if (ii_MinN > n) ii_MinN = n;
            if (ii_MaxN < k) ii_MaxN = n;
         }
      }
   }
   else
      FatalError("Input file %s has unknown format", fileName);
   
   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!memcmp(buffer, "ABCD ", 5))
      {
         if (sscanf(buffer, "ABCD $a*2^%u+1 [%" SCNu64"]", &n, &k) != 2)
            FatalError("Line 1 is not a valid ABCD line in input file %s", fileName);
         
         if (haveBitMap)
         {
            iv_Terms[n-ii_MinN][BIT(k)] = true;
            il_TermCount++;
         }
         else
         {
            if (il_MinK > k) il_MinK = k;
            if (il_MaxK < k) il_MaxK = k;
            if (ii_MinN > n) ii_MinN = n;
            if (ii_MaxN < n) ii_MaxN = n;
         }
         
         continue;
      }

      if (sscanf(buffer, "%" SCNu64 , &diff) != 1)
         FatalError("Line %s is malformed", buffer);
      
      k += diff;
            
      if (haveBitMap)
      {
         iv_Terms[n-ii_MinN][BIT(k)] = true;
         il_TermCount++;
      }
      else
      {
         if (il_MinK > k) il_MinK = k;
         if (il_MaxK < k) il_MaxK = k;
      }
   }
}

bool  GFNDivisorApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint64_t k;
   uint32_t b, n;
   int32_t  c;
   
   if (!ib_UseTermsBitmap)
      return false;
   
   if (sscanf(term, "%" SCNu64"*%u^%u%u", &k, &b, &n, &c) != 4)
      FatalError("Could not parse term %s", term);

   if (b != 2)
      FatalError("Expected base 2 in factor but found base %u", b);
   
   if (c != 1)
      FatalError("Expected c = +1 in factor but found %d", c);
      
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   if (k < il_MinK || k > il_MaxK)
      return false;

   if (!VerifyFactor(false, thePrime, k, n))
      return false;
      
   uint64_t bit = BIT(k);
   
   // No locking is needed because the Workers aren't running yet
   if (iv_Terms[n-ii_MinN][bit])
   {
      iv_Terms[n-ii_MinN][bit] = false;
      il_TermCount--;

      return true;
   }
      
   return false;
}

void GFNDivisorApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   uint32_t fileCount, n;
   char     fileName[200];
   char     termCountStr[50];
   char     termsCountedStr[50];

   if (!ib_UseTermsBitmap)
      return;
      
   // With super large ranges, wait until we can lock because without locking
   // the term count can change between opening and closing the file.
   if (IsRunning() && largestPrime < GetMaxPrimeForSingleWorker())
      return;

   ip_FactorAppLock->Lock();
   
   if (ii_NsPerFile >= (ii_MaxN - ii_MinN + 1))
   {
      sprintf(fileName, "%s.pfgw", is_OutputTermsFilePrefix.c_str());
           
      termsCounted = WriteABCDTermsFile(fileName, ii_MinN, largestPrime);
      
      fileCount = 1;
      
      is_OutputTermsFileName = fileName;
   }
   else
   {   
      fileCount = 0;
      termsCounted = 0;
      
      for (n=ii_MinN; n<=ii_MaxN; n+=ii_NsPerFile)
      {
         fileCount++;
         
         sprintf(fileName, "%s_%04d.pfgw", is_OutputTermsFilePrefix.c_str(), fileCount);
         
         termsCounted += WriteABCDTermsFile(fileName, n, largestPrime);
      }
   }

   if (termsCounted != il_TermCount)
   {
      sprintf(termCountStr, "%" PRIu64"", il_TermCount);
      sprintf(termsCountedStr, "%" PRIu64"", termsCounted);
      
      FatalError("Something is wrong.  Counted terms (%s) != expected terms (%s)", termsCountedStr, termCountStr);
   }
   
   
   ip_FactorAppLock->Release();
}

uint64_t GFNDivisorApp::WriteABCDTermsFile(char *fileName, uint32_t minN, uint64_t maxPrime)
{
   FILE    *termsFile;
   uint32_t n, maxN;
   uint64_t k, bit, termCount = 0, previousK;
   bool     firstRowInFile = true, addedNToFile;

   termsFile = fopen(fileName, "w");
   
   if (!termsFile)
      FatalError("Could not open file %s", fileName);
   
   maxN = minN + ii_NsPerFile;
   if (maxN > ii_MaxN)
      maxN = ii_MaxN + 1;
   
   for (n=minN; n<maxN; n++)
   {
      addedNToFile = false;
      
      k = il_MinK;
      bit = BIT(il_MinK);
   
      for (; k<=il_MaxK; k+=2, bit++)
      {
         if (iv_Terms[n-ii_MinN][bit])
         {
            if (firstRowInFile)
               fprintf(termsFile, "ABCD $a*2^%d+1 [%" PRIu64"] // Sieved to %" PRIu64"\n", n, k, maxPrime);
            else
               fprintf(termsFile, "ABCD $a*2^%d+1 [%" PRIu64"]\n", n, k);

            firstRowInFile = false;
            addedNToFile = true;
            previousK = k;
            termCount++;
            break;
         }
      }

      // If all k for this n were sieved out, we'll skip to the next n.
      if (!addedNToFile)
         continue;
      
      k += 2;
      bit++;
      
      for (; k<=il_MaxK; k+=2, bit++)
      {
         if (iv_Terms[n-ii_MinN][bit])
         {
            fprintf(termsFile, "%" PRIu64"\n", k - previousK);
            previousK = k;
            termCount++;
         }
      }
   }

   fclose(termsFile);
   return termCount;
}

void  GFNDivisorApp::GetExtraTextForSieveStartedMessage(char *extraText)
{
   char  minK[30];
   char  maxK[30];
      
   ConvertNumberToShortString(il_MinK, (char *) minK);
   ConvertNumberToShortString(il_MaxK, (char *) maxK);
   
   sprintf(extraText, "%s <= k <= %s, %d <= n <= %d, k*2^n+1", minK, maxK, ii_MinN, ii_MaxN);
}

bool  GFNDivisorApp::ReportFactor(uint64_t p, uint64_t k, uint32_t n, bool verifyFactor)
{
   bool removedTerm = false;
   char kStr[50];
   char pStr[50];
   
   if (!ib_UseTermsBitmap)
   {
      uint32_t smallPrimeFactor = GetSmallPrimeFactor(k, n);
      
      if (smallPrimeFactor > 0 && smallPrimeFactor <= ii_SmallPrimeFactorLimit)
         return true;

      sprintf(kStr, "%" PRIu64"", k);
      
      il_FactorCount++;
      LogFactor(p, "%s*2^%u+1", kStr, n);

      if (verifyFactor)
         VerifyFactor(true, p, k, n);
         
      return true;
   }
      
   if (p > GetMaxPrimeForSingleWorker() && GetTotalWorkers() > 1)
      ip_FactorAppLock->Lock();

   uint64_t bit = BIT(k);
      
   if (iv_Terms[n-ii_MinN][bit])
   {
      iv_Terms[n-ii_MinN][bit] = false;
      
      il_FactorCount++;
      il_TermCount--;

      sprintf(kStr, "%" PRIu64"", k);
               
      if (n < 62)
      {
         uint64_t nexp = (1L << n);
         
         if (nexp < p)
         {
            if (((p - 1) >> n) == k)
            {
               sprintf(pStr, "%" PRIu64"", p);
               
               WriteToConsole(COT_OTHER, "%s*2^%u+1 is prime (= %s)", kStr, n, pStr);
            }
         }
      }

      LogFactor(p, "%s*2^%u+1", kStr, n);
      removedTerm = true;

      if (verifyFactor)
         VerifyFactor(true, p, k, n);
   }
   
   if (p > GetMaxPrimeForSingleWorker() && GetTotalWorkers() > 1)
      ip_FactorAppLock->Release();

   return removedTerm;
}

uint32_t GFNDivisorApp::GetSmallPrimeFactor(uint64_t k, uint32_t n)
{
   uint32_t idx;
   uint64_t ps[4];
   vector<uint64_t>::iterator it = iv_SmallPrimes.begin();

   while (it != iv_SmallPrimes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;      

      MpArithVec mp(ps);

      const MpResVec pOne = mp.one();
      const MpResVec pTwo = mp.add(mp.one(), mp.one());
      const MpResVec mOne = mp.sub(mp.zero(), pOne);

      MpResVec rem = mp.pow(pTwo, n);  
      rem = mp.mul(rem, mp.nToRes(k));

      for (idx=0; idx<VECTOR_SIZE; idx++)
      {
         if (rem[idx] == mOne[idx])
            return ps[idx];
      }
   }
        
   return 0;
}


bool  GFNDivisorApp::VerifyFactor(bool badFactorIsFatal, uint64_t thePrime, uint64_t k, uint32_t n)
{
   uint64_t rem;

   fpu_push_1divp(thePrime);
   
   rem = fpu_powmod(2, n, thePrime);

   rem = fpu_mulmod(rem, k, thePrime);
   
   fpu_pop();
   
   if (rem == thePrime - 1)
      return true;
      
   char buffer[200];
   
   sprintf(buffer, "Invalid factor: %" PRIu64"*2^%u+1 mod %" PRIu64" = %" PRIu64"", k, n, thePrime, rem+1);
   
   if (badFactorIsFatal)
      FatalError(buffer);
   else
      WriteToConsole(COT_OTHER, buffer);
   
   return false;
}

