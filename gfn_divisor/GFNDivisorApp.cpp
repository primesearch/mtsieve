/* GFNDivisorApp.cpp -- (C) Mark Rodenkirch, November 2017

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
#include "GFNDivisorApp.h"
#include "GFNDivisorWorker.h"
#include "../x86_asm/fpu-asm-x86.h"
#include "../x86_asm_ext/asm-ext-x86.h"

#ifdef HAVE_GPU_WORKERS
#include "GFNDivisorGpuWorker.h"
#endif

#define KMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

#define APP_NAME        "gfndsieve"
#define APP_VERSION     "1.9"

#define BIT(k)          (((k) - il_MinK) >> 1)

// Set PRE_SQUARE=N to compute 2^2^n as (2^2^N)^2^(n-N), which saves N
// sqrmods at a cost of more time in mpn_tdiv_qr().  N must satisfy 0 <= N <= 5.
#define PRE_SQUARE 5

// Handle the possibility that mp_limb_t and uint64_t are different typedefs.
// If they are they the compiler will not be able to cast.
#ifdef WIN32
static const unsigned long long _MONTGOMERY_DATA[14] = {1,0,0,0,0,0,0,0,0,0,0,0,0,UINT64_C(1)<<(1<<PRE_SQUARE)};
#else
static const unsigned long _MONTGOMERY_DATA[14] = {1,0,0,0,0,0,0,0,0,0,0,0,0,UINT64_C(1)<<(1<<PRE_SQUARE)};
#endif

#define ONE    (_MONTGOMERY_DATA+0)
#define TWO128 (_MONTGOMERY_DATA+11)
#define TWO192 (_MONTGOMERY_DATA+10)
#define TWO256 (_MONTGOMERY_DATA+9)
#define TWO320 (_MONTGOMERY_DATA+8)
#define TWO384 (_MONTGOMERY_DATA+7)
#define TWO448 (_MONTGOMERY_DATA+6)
#define TWO512 (_MONTGOMERY_DATA+5)
#define TWO576 (_MONTGOMERY_DATA+4)
#define TWO640 (_MONTGOMERY_DATA+3)
#define TWO704 (_MONTGOMERY_DATA+2)
#define TWO768 (_MONTGOMERY_DATA+1)

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
   il_TotalTerms = 0;
   il_TotalTermsEvaluated = 0;

   SetAppMinPrime(3);
   
#ifdef HAVE_GPU_WORKERS
   ii_GpuFactorDensity = 100;
#endif
}

void GFNDivisorApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k               minimum k to search\n");
   printf("-K --kmax=K               maximum k to search\n");
   printf("-n --nmin=N               minimum n to search\n");
   printf("-N --nmax=N               maximum n to search\n");
   printf("-T --nsperfile=T          number of n per output file\n");
   printf("-x --testterms            test remaining terms for GFN divisibility\n");
   printf("-X --termsperchunk=X      when using -x, number of terms to sieve at a time (default 1e10)\n");
   
#ifdef HAVE_GPU_WORKERS
   printf("-M --maxfactordensity=M   factors per 1e6 terms per GPU worker chunk (default %u)\n", ii_GpuFactorDensity);
#endif
}

void  GFNDivisorApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "k:K:n:N:c:T:xX:";

   AppendLongOpt(longOpts, "kmin",              required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",              required_argument, 0, 'K');
   AppendLongOpt(longOpts, "nmin",              required_argument, 0, 'n');
   AppendLongOpt(longOpts, "nmax",              required_argument, 0, 'N');
   AppendLongOpt(longOpts, "nsperfile",         required_argument, 0, 'T');
   AppendLongOpt(longOpts, "testterms",         no_argument,       0, 'x');
   AppendLongOpt(longOpts, "termsperchunk",     required_argument, 0, 'X');
   
#ifdef HAVE_GPU_WORKERS
   shortOpts += "M:";
   
   AppendLongOpt(longOpts, "maxfactordensity",  required_argument, 0, 'M');
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
      case 'M':
         status = Parser::Parse(arg, 10, 10000000, ii_GpuFactorDensity);
         break;
#endif
   }

   return status;
}

void GFNDivisorApp::ValidateOptions(void)
{
   uint64_t kCount;
   uint32_t nCount;
   
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
      
      iv_Terms.resize(nCount);
      
      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         iv_Terms[n-ii_MinN].resize(kCount);
         std::fill(iv_Terms[n-ii_MinN].begin(), iv_Terms[n-ii_MinN].end(), false);
      }
            
      ProcessInputTermsFile(true);
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

   // Not expecting anyone to use gfndsieve for n > 2000 since gfndsieve+pfgw should be faster
   if (ib_TestTerms && il_MaxPrime == il_AppMaxPrime) {
      if (ii_MaxN <= 100000) SetMaxPrime(30000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=  50000) SetMaxPrime(30000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=  30000) SetMaxPrime(30000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=  20000) SetMaxPrime(30000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=  10000) SetMaxPrime(30000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=   5000) SetMaxPrime(30000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=   3000) SetMaxPrime(10000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=   2000) SetMaxPrime(10000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=   1000) SetMaxPrime(10000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    800) SetMaxPrime(10000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    600) SetMaxPrime(10000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    400) SetMaxPrime( 3000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    300) SetMaxPrime( 3000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    200) SetMaxPrime( 1000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    150) SetMaxPrime( 1000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=    100) SetMaxPrime( 1000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=     80) SetMaxPrime( 1000000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=     60) SetMaxPrime(  100000000, "testing terms for GFN divisibility");
      if (ii_MaxN <=     40) SetMaxPrime(    3000000, "testing terms for GFN divisibility");
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
   
#ifdef HAVE_GPU_WORKERS
   double factors = (double) (ii_MaxN - ii_MinN) * (double) (il_MaxK - il_MinK) / 1000000.0;

   ii_MaxGpuFactors = GetGpuWorkGroups() * (uint64_t) (factors * (double) ii_GpuFactorDensity);
#endif
}

void  GFNDivisorApp::PreSieveHook(void)
{
   uint64_t kCount;
   uint32_t nCount;
   
   kCount = (il_MaxK - il_MinK)/2 + 1;
   nCount = ii_MaxN - ii_MinN + 1;
   il_TotalTerms = kCount * nCount;
   
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
   
   il_StartSievingUS = Clock::GetCurrentMicrosecond();
}

bool  GFNDivisorApp::PostSieveHook(void)
{
   if (!ib_TestTerms)
      return true;
   
   TestRemainingTerms();
   
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
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);
   
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
      FatalError("Could not open file %s", termsFile);
   
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
   if (k < il_MinK || k > il_MaxK)
      return false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   bool removedTerm = false;
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Lock();

   uint64_t bit = BIT(k);
      
   if (iv_Terms[n-ii_MinN][bit])
   {
      iv_Terms[n-ii_MinN][bit] = false;
      
      il_FactorCount++;
      il_TermCount--;

      if (n < 62)
      {
         uint64_t nexp = (1L << n);
         
         if (nexp < p)
         {
            if (((p - 1) >> n) == k)
               WriteToConsole(COT_OTHER, "%" PRIu64"*2^%u+1 is prime (= %" PRIu64")", k, n, p);
         }
      }

      LogFactor(p, "%" PRIu64"*2^%u+1", k, n);
      removedTerm = true;

      if (verifyFactor)
         VerifyFactor(true, p, k, n);
   }
   
   if (p > GetMaxPrimeForSingleWorker())
      ip_FactorAppLock->Release();

   return removedTerm;
}

void  GFNDivisorApp::TestRemainingTerms(void)
{
   time_t   lastCheckPointTime = time(NULL);
   uint64_t sievingUS, startTestingUS, currentUS;
   uint64_t calculationUS;
   uint64_t termsTested = 0;
   uint64_t termsEvaluatedInChunk = 0;
   uint64_t termsTestedPerSecond;
   double   percentSievingTimeSlice;
   double   percentTermsRequiringTest;
   double   percentChunkTested, percentRangeTested;
   FILE    *fPtr;
   mpz_t    rem, fermat, nTemp, kTemp, factor, minus1;

   startTestingUS = Clock::GetCurrentMicrosecond();
   sievingUS = startTestingUS - il_StartSievingUS;
            
   mpz_init(rem);
   mpz_init(fermat);
   mpz_init(nTemp);
   mpz_init(kTemp);
   mpz_init(factor);
   mpz_init(minus1);

   mpz_set_ui(fermat, 2);
   
   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      uint64_t bit = BIT(il_MinK);
         
      mpz_set_ui(nTemp, 2);
      mpz_pow_ui(nTemp, nTemp, n);
         
      for (uint64_t k=il_MinK; k<=il_MaxK; k+=2, bit++)
      {
         il_TotalTermsEvaluated++;
         termsEvaluatedInChunk++;
         
         if (!iv_Terms[n-ii_MinN][bit])
            continue;

         termsTested++;

         if (time(NULL) > lastCheckPointTime + 60)
         {            
            percentSievingTimeSlice = ((double) termsEvaluatedInChunk) / (double) il_TotalTermsInChunk;
            percentChunkTested = 100.0 * percentSievingTimeSlice;
            
            currentUS = Clock::GetCurrentMicrosecond();
            
            // Add the time to test this range to the time to sieve this range.  So if it took 180 seconds
            // to sieve this chunk and we have tested 20 percent of this chunk then "assign" 36 seconds
            // (as 36 is 20% of 180) of sieving time to this chunk.
            calculationUS = (sievingUS * percentSievingTimeSlice) + (currentUS - startTestingUS);
                        
            percentTermsRequiringTest = (100.0 * (double) il_TermCount) / (double) il_TotalTermsInChunk;
            
            termsTestedPerSecond = termsEvaluatedInChunk / (calculationUS / 1000000);
            
            // We really didn't evaluate even k, but we count against the rate anyways.
            termsTestedPerSecond *= 2;

            if (il_TotalTermsInChunk == il_TotalTerms)
               WriteToConsole(COT_SIEVE, "Tested %5.2f pct of range at %" PRIu64" terms per second (%5.2f pct terms passed sieving)", percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest);
            else
            {
               percentRangeTested = (100.0 * (double) il_TotalTermsEvaluated) / (double) il_TotalTerms;
            
               WriteToConsole(COT_SIEVE, "Tested %5.2f pct of chunk at %" PRIu64" terms per second (%5.2f pct terms passed sieving) (%5.2f pct of range)", percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest, percentRangeTested);
            }
      
            lastCheckPointTime = time(NULL);
         }
         
         if (IsFermatDivisor(k, n))
         {
#ifdef WIN32
            // Even though build with 64-bit limbs, mpz_set_ui doesn't
            // populate kTemp correctly when k > 32 bits.
            mpz_set_ui(kTemp, k >> 32);
            mpz_mul_2exp(kTemp, kTemp, 32);
            mpz_add_ui(kTemp, kTemp, k & (0xffffffff));
#else
            mpz_set_ui(kTemp, k);
#endif

            mpz_mul_2exp(minus1, kTemp, n);
            mpz_add_ui(factor, minus1, 1);

            mpz_powm(rem, fermat, nTemp, factor);;
         
            if (mpz_cmp(rem, minus1) == 0)
            {
               WriteToConsole(COT_OTHER, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1", k, n, n-2);
               
               fPtr = fopen("gfn_factors.txt", "a+");
               fprintf(fPtr, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1\n", k, n, n-2);
               fclose(fPtr);
            }
            else if (mpz_cmp_ui(rem, 1) == 0)
            {
               mpz_set_ui(rem, 2);
               
               for (uint32_t m=1; m<=n-2; m++)
               {
                  mpz_powm_ui(rem, rem, 2, factor);
                  
                  if (mpz_cmp(rem, minus1) == 0)
                  {
                     WriteToConsole(COT_OTHER, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1", k, n, m);
                     
                     fPtr = fopen("gfn_factors.txt", "a+");
                     fprintf(fPtr, "Found factor %" PRIu64"*2^%u+1 of 2^(2^%u)+1\n", k, n, m);
                     fclose(fPtr);
                  }
               }
            }
         }
      }
   }
   
   mpz_clear(rem);
   mpz_clear(fermat);
   mpz_clear(nTemp);
   mpz_clear(kTemp);
   mpz_clear(factor);
   mpz_clear(minus1);

   currentUS = Clock::GetCurrentMicrosecond();
   
   // If we took less than 5 seconds to test the range, then the terms per second
   // calculation is rather meaningless so we won't show it.
   if (currentUS - il_StartSievingUS > 50000000)
   {
      percentChunkTested = (100.0 * (double) termsEvaluatedInChunk) / (double) il_TotalTermsInChunk;
      percentTermsRequiringTest = (100.0 * (double) il_TermCount) / (double) il_TotalTermsInChunk;
      
      termsTestedPerSecond = termsEvaluatedInChunk / ((currentUS - il_StartSievingUS) / 1000000);

      // We really didn't evaluate even k, but we count against the rate anyways.
      termsTestedPerSecond *= 2;
               
      if (il_TotalTermsInChunk == il_TotalTerms)
         WriteToConsole(COT_SIEVE, "Tested %5.2f pct of range at %" PRIu64" terms per second (%5.2f pct terms passed sieving)", percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest);
      else
      {
         percentRangeTested = (100.0 * (double) il_TotalTermsEvaluated) / (double) il_TotalTerms;
      
         WriteToConsole(COT_SIEVE, "Tested %5.2f pct of chunk at %" PRIu64" terms per second (%5.2f pct terms passed sieving) (%5.2f pct of range)", percentChunkTested, termsTestedPerSecond, percentTermsRequiringTest, percentRangeTested);
      }
   }
}

// Return 1 iff k*2^n+1 is a Fermat divisor.
// This code is from fermat_redc.c of GMP-Fermat
bool  GFNDivisorApp::IsFermatDivisor(uint64_t k, uint32_t n)
{
// Handle the possibility that mp_limb_t and uint64_t are different typedefs.
// We don't use mp_limb_t here, but the variables are passed to our assembler
// ext functions and to GMP functions, so theu must be compatible.
#ifdef WIN32
   unsigned long long inv, A[14], N[12];
#else
   unsigned long inv, A[14], N[12];
#endif

   uint32_t i;


   // If 2^2^(n-1) = 1 (mod k*2^n+1) then k*2^n+1 divides 2^2^(n-1)-1 and
   // hence divides at least one of 2^2^1+1, 2^2^2+1, ..., 2^2^(n-2)+1.
   switch (n/64)
   {
      case 0: /* n < 64, so k*2^n+1 < 2^127 */
         /* N <-- k*2^n+1 */
         N[0] = (k << n) | 1;
         N[1] = (k >> (64-n));

         /* Encode 2 into Montgomery form, A <-- 2*2^128 mod N */
         A[1] = 0;
         
         mpn_tdiv_qr(A+2,A,0L,TWO128,3,N,(N[1]>0)?2:1);

         /* A <-- 2^2^(n-1) mod N, in Montgomery form */
         inv = invproth(k,n);
         for (i = n-1-PRE_SQUARE; i > 0; i--)
            sqrmod128_proth0(A,N,inv);

         /* Decode result from Montgomery form */
         mulmod128(A,ONE,A,N,inv);
         CheckRedc(A,2,2,n-1,k,n);

         /* Return 1 iff k*2^n+1 is a Fermat divisor */
         return (A[0] == 1 && A[1] == 0);

      case 1: /* n < 128 */
         if (n == 64 || (k >> (128-n)) == 0) /* 2^64 < k*2^n+1 < 2^128 */
         {
            N[0] = 1;
            N[1] = (k << (n-64));
            
            mpn_tdiv_qr(A+2,A,0L,TWO128,3,N,2);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod128_proth1(A,N);
         
            mulmod128_proth1(A,ONE,A,N);
            CheckRedc(A,2,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0);
         }
         else /* 2^128 < k*2^n+1 < 2^191 */
         {
            N[0] = 1;
            N[1] = (k << (n-64));
            N[2] = (k >> (128-n));
            
            mpn_tdiv_qr(A+3,A,0L,TWO192,4,N,3);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod192_proth0(A,N);
         
            mulmod192_proth0(A,ONE,A,N);
            
            CheckRedc(A,3,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0);
         }

      case 2: /* n < 192 */
         if (n == 128 || (k >> (192-n)) == 0) /* 2^128 < k*2^n+1 < 2^192 */
         {
            N[0] = 1; N[1] = 0;
            N[2] = (k << (n-128));
            
            mpn_tdiv_qr(A+3,A,0L,TWO192,4,N,3);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod192_proth1(A,N);
         
            mulmod192_proth1(A,ONE,A,N);
            
            CheckRedc(A,3,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0);
         }
         else /* 2^192 < k*2^n+1 < 2^255 */
         {
            N[0] = 1; N[1] = 0;
            N[2] = (k << (n-128));
            N[3] = (k >> (192-n));
            
            mpn_tdiv_qr(A+4,A,0L,TWO256,5,N,4);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod256_proth0(A,N);
            
            mulmod256_proth0(A,ONE,A,N);
            
            CheckRedc(A,4,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0);
         }

      case 3: /* n < 256 */
         if (n == 192 || (k >> (256-n)) == 0) /* 2^192 < k*2^n+1 < 2^256 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0;
            N[3] = (k << (n-192));
            
            mpn_tdiv_qr(A+4,A,0L,TWO256,5,N,4);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               sqrmod256_proth1(A,N);
            
            mulmod256_proth1(A,ONE,A,N);
            
            CheckRedc(A,4,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0);
         }
         else /* 2^256 < k*2^n+1 < 2^319 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0;
            N[3] = (k << (n-192));
            N[4] = (k >> (256-n));
            
            mpn_tdiv_qr(A+5,A,0L,TWO320,6,N,5);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod320_proth0(A,A,A,N);
         
            mulmod320_proth0(A,ONE,A,N);
            
            CheckRedc(A,5,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0);
         }

      case 4: /* n < 320 */
         if (n == 256 || (k >> (320-n)) == 0) /* 2^256 < k*2^n+1 < 2^320 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0;
            N[4] = (k << (n-256));
            
            mpn_tdiv_qr(A+5,A,0L,TWO320,6,N,5);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod320_proth1(A,A,A,N);
         
            mulmod320_proth1(A,ONE,A,N);
            
            CheckRedc(A,5,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0);
         }
         else /* 2^320 < k*2^n+1 < 2^383 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0;
            N[4] = (k << (n-256));
            N[5] = (k >> (320-n));
            
            mpn_tdiv_qr(A+6,A,0L,TWO384,7,N,6);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod384_proth0(A,A,A,N);
            
            mulmod384_proth0(A,ONE,A,N);
            
            CheckRedc(A,6,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0);
         }

      case 5: /* n < 384 */
         if (n == 320 || (k >> (384-n)) == 0) /* 2^320 < k*2^n+1 < 2^384 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0;
            N[5] = (k << (n-320));
            
            mpn_tdiv_qr(A+6,A,0L,TWO384,7,N,6);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod384_proth1(A,A,A,N);
         
            mulmod384_proth1(A,ONE,A,N);
            
            CheckRedc(A,6,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0);
         }
         else /* 2^384 < k*2^n+1 < 2^447 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0;
            N[5] = (k << (n-320));
            N[6] = (k >> (384-n));
            
            mpn_tdiv_qr(A+7,A,0L,TWO448,8,N,7);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod448_proth0(A,A,A,N);
         
            mulmod448_proth0(A,ONE,A,N);
            
            CheckRedc(A,7,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0 && A[6] == 0);
         }

      case 6: /* n < 448 */
         if (n == 384 || (k >> (448-n)) == 0) /* 2^384 < k*2^n+1 < 2^448 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = (k << (n-384));
            
            mpn_tdiv_qr(A+7,A,0L,TWO448,8,N,7);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod448_proth1(A,A,A,N);
         
            mulmod448_proth1(A,ONE,A,N);
            
            CheckRedc(A,7,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && A[5] == 0 && A[6] == 0);
         }
         else /* 2^448 < k*2^n+1 < 2^511 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = (k << (n-384));
            N[7] = (k >> (448-n));
            
            mpn_tdiv_qr(A+8,A,0L,TWO512,9,N,8);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod512_proth0(A,A,A,N);
         
            mulmod512_proth0(A,ONE,A,N);
            
            CheckRedc(A,8,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0);
         }

      case 7: /* n < 512 */
         if (n == 448 || (k >> (512-n)) == 0) /* 2^448 < k*2^n+1 < 2^512 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0;
            N[7] = (k << (n-448));
            
            mpn_tdiv_qr(A+8,A,0L,TWO512,9,N,8);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod512_proth1(A,A,A,N);
         
            mulmod512_proth1(A,ONE,A,N);
            
            CheckRedc(A,8,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0);
         }
         else /* 2^512 < k*2^n+1 < 2^575 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0;
            N[7] = (k << (n-448));
            N[8] = (k >> (512-n));
            
            mpn_tdiv_qr(A+9,A,0L,TWO576,10,N,9);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod576_proth0(A,A,A,N);
         
            mulmod576_proth0(A,ONE,A,N);
            
            CheckRedc(A,9,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0);
         }

      case 8: /* n < 576 */
         if (n == 512 || (k >> (576-n)) == 0) /* 2^512 < k*2^n+1 < 2^576 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0;
            N[8] = (k << (n-512));
            
            mpn_tdiv_qr(A+9,A,0L,TWO576,10,N,9);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod576_proth1(A,A,A,N);
         
            mulmod576_proth1(A,ONE,A,N);
            
            CheckRedc(A,9,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0);
         }
         else /* 2^576 < k*2^n+1 < 2^639 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0;
            N[8] = (k << (n-512));
            N[9] = (k >> (576-n));
            
            mpn_tdiv_qr(A+10,A,0L,TWO640,11,N,10);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod640_proth0(A,A,A,N);
         
            mulmod640_proth0(A,ONE,A,N);
            
            CheckRedc(A,10,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 && A[9] == 0);
         }

      case 9: /* n < 640 */

         if (n == 576 || (k >> (640-n)) == 0) /* 2^576 < k*2^n+1 < 2^640 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0;
            N[9] = (k << (n-576));

            mpn_tdiv_qr(A+10,A,0L,TWO640,11,N,10);

            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod640_proth1(A,A,A,N);

            mulmod640_proth1(A,ONE,A,N);
            
            CheckRedc(A,10,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 && 
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 && A[9] == 0);
         }
         else /* 2^640 < k*2^n+1 < 2^703 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0;
            N[9] = (k << (n-576));
            N[10] = (k >> (640-n));
            
            mpn_tdiv_qr(A+11,A,0L,TWO704,12,N,11);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod704_proth0(A,A,A,N);
         
            mulmod704_proth0(A,ONE,A,N);
            
            CheckRedc(A,11,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 && A[9] == 0 && A[10] == 0);
         }

      case 10: /* n < 704 */
         if (n == 640 || (k >> (704-n)) == 0) /* 2^640 < k*2^n+1 < 2^704 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0; N[9] = 0;
            N[10] = (k << (n-640));
            
            mpn_tdiv_qr(A+11,A,0L,TWO704,12,N,11);

            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod704_proth1(A,A,A,N);

            mulmod704_proth1(A,ONE,A,N);

            CheckRedc(A,11,2,n-1,k,n);

            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 &&
                    A[9] == 0 && A[10] == 0);
         }
         else /* 2^704 < k*2^n+1 < 2^767 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0;
            N[6] = 0; N[7] = 0; N[8] = 0; N[9] = 0;
            N[10] = (k << (n-640));
            N[11] = (k >> (704-n));

            mpn_tdiv_qr(A+12,A,0L,TWO768,13,N,12);

            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod768_proth0(A,A,A,N);

            mulmod768_proth0(A,ONE,A,N);

            CheckRedc(A,12,2,n-1,k,n);

            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 &&
                    A[9] == 0 && A[10] == 0 && A[11] == 0);
         }

      case 11: /* n < 768 */
         if (n == 704 || (k >> (768-n)) == 0) /* 2^704 < k*2^n+1 < 2^768 */
         {
            N[0] = 1; N[1] = 0; N[2] = 0; N[3] = 0; N[4] = 0; N[5] = 0; N[6] = 0;
            N[7] = 0; N[8] = 0; N[9] = 0; N[10] = 0;
            N[11] = (k << (n-704));
            
            mpn_tdiv_qr(A+12,A,0L,TWO768,13,N,12);
            
            for (i = n-1-PRE_SQUARE; i > 0; i--)
               mulmod768_proth1(A,A,A,N);
         
            mulmod768_proth1(A,ONE,A,N);
            
            CheckRedc(A,12,2,n-1,k,n);
            
            return (A[0] == 1 && A[1] == 0 && A[2] == 0 && A[3] == 0 && A[4] == 0 &&
                    A[5] == 0 && A[6] == 0 && A[7] == 0 && A[8] == 0 &&
                    A[9] == 0 && A[10] == 0 && A[11] == 0);
         }

      default:
      {
         uint32_t nn = n/64;
         mp_limb_t np[nn+2], t1[2*nn+4], t2[2*nn+4];

         np[0] = 1;
         
         for (i = 1; i < nn; i++)
            np[i] = 0;
         
         np[nn++] = k<<(n%64);
         
         if (n%64)
            if ((np[nn] = k>>(64-n%64)) != 0)
               nn++;

         for (i = 0; i < nn; i++)
            t1[i] = 0;
         
         t1[nn] = (mp_limb_t)1<<(1<<PRE_SQUARE);
         mpn_tdiv_qr(t2+nn,t2,0L,t1,nn+1,np,nn);

         for (i = n-1-PRE_SQUARE; i > 0; i--)
         {
            mpn_mul_n(t1,t2,t2,nn);
            redc_proth0(t1,np,nn);
            if (mpn_add_n(t2,t1+nn,t1,nn))
               mpn_sub_n(t2,t2,np,nn);
         }
         
         for (i = nn; i < 2*nn; i++)
            t2[i] = 0;
         
         redc_proth0(t2,np,nn);
         
         if (mpn_add_n(t1,t2+nn,t2,nn))
            mpn_sub_n(t1,t1,np,nn);
      
#if 1
         // Is this needed?
         while(mpn_cmp(t1,np,nn) >= 0)
            mpn_sub_n(t1,t1,np,nn);
#endif
         
         CheckRedc(t1,nn,2,n-1,k,n);
         
         if (t1[0] != 1)
            return false;
      
         for (i = 1; i < nn; i++)
            if (t1[i] != 0)
               return false;
         
         return true;
      }
   }
}

void  GFNDivisorApp::CheckRedc(mp_limb_t *xp, uint32_t xn, uint32_t b, uint32_t m, uint64_t k, uint32_t n)
{
#ifndef TEST_REDC
   return;
#endif

   mpz_t     B, E, N;
   uint32_t  i;

   mpz_init_set_ui(B,b);

   mpz_init_set_ui(E,1);
   mpz_mul_2exp(E,E,m);

   mpz_init(N);
   
#ifdef WIN32
   // Even though build with 64-bit limbs, mpz_set_ui doesn't
   // populate kTemp correctly when k > 32 bits.
   mpz_set_ui(N, k >> 32);
   mpz_mul_2exp(N, N, 32);
   mpz_add_ui(N, N, k & (0xffffffff));
#else
   mpz_set_ui(N, k);
#endif

   mpz_mul_2exp(N,N,n);
   mpz_add_ui(N,N,1);

   mpz_powm(B,B,E,N);

   if (mpz_size(B) > xn)
   {
      WriteToConsole(COT_OTHER, "REDC ERROR: %d^2^%d mod %" PRIu64"*2^%d+1\n",b,m,k,n);
      FatalError("mpz_size=%lu,xn=%lu\n",mpz_size(B),xn);
   }

   for (i = 0; i < xn; i++)
      if (xp[i] != mpz_getlimbn(B,i))
      {
         WriteToConsole(COT_OTHER, "REDC ERROR: %d^2^%d mod %" PRIu64"*2^%d+1\n",b,m,k,n);
         FatalError("limb %d=%lu, xp[%d]=%lu\n",i,mpz_getlimbn(B,i),i,xp[i]);
      }

   mpz_clear(N);
   mpz_clear(E);
   mpz_clear(B);
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

