/* TwinApp.cpp -- (C) Mark Rodenkirch, September 2018

   Sieve for k*b^n+1 and k*b^n-1 for a range of k and a fixed b and n.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include "../core/Parser.h"
#include "../core/Clock.h"
#include "../core/MpArith.h"
#include "TwinApp.h"
#include "TwinWorker.h"

#define APP_NAME        "twinsieve"
#define APP_VERSION     "1.6"

#define NMAX_MAX        (1 << 31)
#define BMAX_MAX        (1 << 31)

#define BIT_HK(k)       (((k) - il_MinK) >> 1)
#define BIT_AK(k)       ((k) - il_MinK)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the mtsieve framework for other applications.
App *get_app(void)
{
   return new TwinApp();
}

TwinApp::TwinApp() : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors of k*b^n+1/-1 numbers for fixed b and n and variable k");
   SetLogFileName("twinsieve.log");
   
   il_MinK = 0;
   il_MaxK = 0;
   ii_Base = 0;
   ii_N    = 0;
   ib_OnlyTwins = true;
   it_Format = FF_ABCD;
   ib_Remove = false;
   it_TermType = TT_UNKNOWN;
   il_Terms = NULL;
   ii_Primes = NULL;
   ib_HalfK = false;
   
   il_MaxPrimeForValidFactor = PMAX_MAX_62BIT;
   
   iv_TwinTerms.clear();
   iv_MinusTerms.clear();
   iv_PlusTerms.clear();
}

TwinApp::~TwinApp()
{
   if (ii_Primes != 0)
      xfree(ii_Primes);
   
   if (il_Terms != 0)
      xfree(il_Terms);
}

void TwinApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-k --kmin=k           Minimum k to search\n");
   printf("-K --kmax=K           Maximum k to search\n");
   printf("-b --base=b           Base to search\n");
   printf("-n --n=n              n of b^n, n# for primorial, n! for factorial\n");
   printf("-t --termtype=t       1 = b^n, 2 = primorial, 3 = factorial\n");
   printf("-f --format=f         Format of output file (A=ABC, D=ABCD (default), N=NEWPGEN)\n");
   printf("-r --remove           Remove k where k %% base = 0\n");
   printf("-s --independent      Sieve +1 and -1 independently\n");       
}

void  TwinApp::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "srk:K:b:n:t:f:";

   AppendLongOpt(longOpts, "kmin",           required_argument, 0, 'k');
   AppendLongOpt(longOpts, "kmax",           required_argument, 0, 'K');
   AppendLongOpt(longOpts, "base",           required_argument, 0, 'b');
   AppendLongOpt(longOpts, "exp",            required_argument, 0, 'n');
   AppendLongOpt(longOpts, "termtype",       required_argument, 0, 't');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
   AppendLongOpt(longOpts, "remove",         no_argument, 0, 'r');
   AppendLongOpt(longOpts, "independent",    no_argument, 0, 's');
}

parse_t TwinApp::ParseOption(int opt, char *arg, const char *source)
{
   char value;
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
         status = Parser::Parse(arg, 2, BMAX_MAX, ii_Base);
         break;
         
      case 'n':
         status = Parser::Parse(arg, 1, NMAX_MAX, ii_N);
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
         
      case 'f':
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
         ib_OnlyTwins = false;
         status = P_SUCCESS;
         break;
         
      case 'r':
         ib_Remove = true;
         status = P_SUCCESS;
         break;
   }

   return status;
}

void TwinApp::ValidateOptions(void)
{
   if (it_Format == FF_UNKNOWN)
      FatalError("File format not valid, use A (ABC), D (ABCD) or N (NewPGen)");

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      // We only care about odd k or even k, but not both
      if (it_TermType == TT_BN && (ii_Base == 2 || ii_Base & 1))
         ib_HalfK = true;
   
      if (it_TermType == TT_PRIMORIAL)
         BuildPrimorialTerms();
      
      if (it_TermType == TT_FACTORIAL)
         BuildFactorialTerms();
      
      uint32_t termSize = (il_MaxK - il_MinK) + 1;
      
      if (ib_HalfK)
         termSize = (il_MaxK - il_MinK)/2 + 1;
            
      if (ib_OnlyTwins)
      {
         iv_TwinTerms.resize(termSize);
         std::fill(iv_TwinTerms.begin(), iv_TwinTerms.end(), false);
      }
      else 
      {
         iv_MinusTerms.resize(termSize);
         std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
         
         iv_PlusTerms.resize(termSize);         
         std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
      }
      
      ProcessInputTermsFile(true);
   }
   else
   {
      if (it_TermType == TT_UNKNOWN)
         FatalError("Term type must be specified");
      
      if (il_MinK == 0)
         FatalError("kmin must be specified");

      if (il_MaxK == 0)
         FatalError("kmax must be specified");
      
      if (il_MaxK <= il_MinK)
         FatalError("kmax must be greater than kmin");
      
      if (it_TermType == TT_BN && ii_Base == 0)
         FatalError("base must be specified");
      
      if (ii_N == 0)
         FatalError("exponent must be specified");

      if (ib_OnlyTwins)
         il_TermCount = (il_MaxK - il_MinK) + 1;
      else
         il_TermCount = 2*((il_MaxK - il_MinK) + 1);
         
      if (it_TermType == TT_BN)
      {
         if (il_MaxK < ii_Base)
            FatalError("base must be less then maxk");
         
         if (ii_Base == 2)
         {
            // We only care about odd k
            ib_HalfK = true;
            if (ib_OnlyTwins)
               il_TermCount = (il_MaxK - il_MinK)/2 + 1;
            else
               il_TermCount = 2*((il_MaxK - il_MinK)/2 + 1);
         
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
            if (ib_OnlyTwins)
               il_TermCount = (il_MaxK - il_MinK)/2 + 1;
            else
               il_TermCount = 2*((il_MaxK - il_MinK)/2 + 1);
         
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
      
      if (ib_OnlyTwins)
      {
         iv_TwinTerms.resize(il_TermCount);
         std::fill(iv_TwinTerms.begin(), iv_TwinTerms.end(), true);
      }
      else 
      {
         iv_MinusTerms.resize(il_TermCount);
         std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), true);
         
         iv_PlusTerms.resize(il_TermCount);
         std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), true);
      } 
   }
   
   if (!ib_OnlyTwins && it_Format != FF_ABC)
   {
      it_Format = FF_ABC;
      WriteToConsole(COT_OTHER, "Switching to ABC format since other formats are not supported when using -s");
   }
   
   
   if (it_TermType != TT_BN && it_Format == FF_NEWPGEN)
   {
      it_Format = FF_ABC;
      WriteToConsole(COT_OTHER, "Switching to ABC format since newpgen format is not supported for primorials or factorials");
   }
         
   
   if (ib_Remove && it_TermType == TT_BN)
   {
      uint64_t k = il_MinK;
      
      while (k % ii_Base > 0)
         k++;
      
      // Remvoe k that are divisible by the base
      for ( ; k<=il_MaxK; k+=ii_Base)
      {
         uint64_t bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));
         
         if (ib_OnlyTwins)
         {
            if (iv_TwinTerms[bit])
            {
               il_TermCount--;
               iv_TwinTerms[bit] = false;
            }
         }
         else
         {
            if (iv_MinusTerms[bit])
            {
               il_TermCount--;
               iv_MinusTerms[bit] = false;
            }
            
            if (iv_PlusTerms[bit])
            {
               il_TermCount--;
               iv_PlusTerms[bit] = false;
            }
         }
      }
   }
         
   if (is_OutputTermsFileName.length() == 0)
   {
      char  fileName[30];
      
      if (it_TermType == TT_BN)
      {
         if (it_Format == FF_NEWPGEN)
            snprintf(fileName, sizeof(fileName), "k_b%u_n%u.npg", ii_Base, ii_N);
         else
            snprintf(fileName, sizeof(fileName), "k_b%u_n%u.pfgw", ii_Base, ii_N);
      }
            
      if (it_TermType == TT_PRIMORIAL)
         snprintf(filePrefix, sizeof(fileName), "twin_%up", ii_N);
      
      if (it_TermType == TT_FACTORIAL)
         snprintf(filePrefix, sizeof(fileName), "twin_%uf", ii_N);
      
      is_OutputTermsFileName = fileName;
   }


   FactorApp::ParentValidateOptions();

   // Since the worker wants primes in groups of 4
   while (ii_CpuWorkSize % 4 != 0)
      ii_CpuWorkSize++;
   
   AdjustMaxPrime();
   
   // Allow only one worker to do work when processing small primes.  This allows us to avoid 
   // locking when factors are reported, which significantly hurts performance as most terms 
   // will be removed due to small primes.
   SetMaxPrimeForSingleWorker(10000);
}

Worker *TwinApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

   // Note that Twin inherits from Worker.  This will not
   // only create the worker, but also start it.
   theWorker = new TwinWorker(id, this);

   return theWorker;
}

void TwinApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   int32_t  c;
   uint32_t n;
   uint32_t base1, base2, n1, n2;
   uint64_t bit, k, diff, lastPrime = 0;
   format_t format = FF_UNKNOWN;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("No data in input file %s", is_InputTermsFileName.c_str());
   
   if (!haveBitMap)
      il_MinK = il_MaxK = 0;
   
   it_TermType = TT_UNKNOWN;
   
   if (!memcmp(buffer, "ABCD ", 5))
   {
      if (strchr(buffer, '&') != NULL)
      {
         if (sscanf(buffer, "ABCD $a*%u^%d+1 & $a*%u^%d-1  [%" SCNu64"] // Sieved to %" SCNu64"", &base1, &n1, &base2, &n2, &k, &lastPrime) == 6)
            it_TermType = TT_BN;

         if (sscanf(buffer, "ABCD $a*%u#+1 & $a*%u#-1  [%" SCNu64"] // Sieved to %" SCNu64"", &n1, &n2, &k, &lastPrime) == 4)
            it_TermType = TT_PRIMORIAL;

         if (sscanf(buffer, "ABCD $a*%u!+1 & $a*%u!-1  [%" SCNu64"] // Sieved to %" SCNu64"", &n1, &n2, &k, &lastPrime) == 4)
            it_TermType = TT_FACTORIAL;

         if (it_TermType == TT_UNKNOWN)
            FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());
            
         if (it_TermType == TT_BN && base1 != base2)
            FatalError("Line 1 bases in ABCD input file %s do not match", is_InputTermsFileName.c_str());
            
         if (n1 != n2)
            FatalError("Line 1 exponents in ABCD input file %s do not match", is_InputTermsFileName.c_str());
         
         ii_Base = base1;
         ii_N = n1;
         ib_OnlyTwins = true;
      }
      else
         FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
            
      format = FF_ABCD;

      if (haveBitMap)
      {
         iv_TwinTerms[k-il_MinK] = true;
         il_TermCount++;
      }
      else
         il_MinK = il_MaxK = k;
   }
   else if (!memcmp(buffer, "ABC ", 4))
   {
      if (strchr(buffer, '&') != NULL)
      {
         if (sscanf(buffer, "ABC $a*%u^%d+1 & $a*%u^%d-1  // Sieved to %" SCNu64"", &base1, &n1, &base2, &n2, &lastPrime) == 5)
            it_TermType = TT_BN;
            
         if (sscanf(buffer, "ABC $a*%u#+1 & $a*%u#-1  // Sieved to %" SCNu64"", &n1, &n2, &lastPrime) == 3)
            it_TermType = TT_PRIMORIAL;

         if (sscanf(buffer, "ABC $a*%u!+1 & $a*%u!-1  // Sieved to %" SCNu64"", &n1, &n2, &lastPrime) == 3)
            it_TermType = TT_FACTORIAL;

         if (it_TermType == TT_UNKNOWN)
            FatalError("Line 1 is not a valid ABCD line in input file %s", is_InputTermsFileName.c_str());
         
         if (it_TermType == TT_BN && base1 != base2)
            FatalError("Line 1 bases in ABC input file %s do not match", is_InputTermsFileName.c_str());
            
         if (n1 != n2)
            FatalError("Line 1 exponents in ABC input file %s do not match", is_InputTermsFileName.c_str());
         
         ii_Base = base1;
         ii_N = n1;
         ib_OnlyTwins = true;
      }
      else
      {
         if (sscanf(buffer, "ABC $a*%u^%d$b // Sieved to %" SCNu64"", &ii_Base, &ii_N, &lastPrime) != 3)
            FatalError("Line 1 is not a valid ABC line in input file %s", is_InputTermsFileName.c_str());
         
         ib_OnlyTwins = false;
      }
      
      format = FF_ABC;
   }
   else if (isdigit(buffer[0]))
   {
      // handle newpgnen noise created by srsieve/srfile
      if (sscanf(buffer, "%" SCNu64":T:0:%u:3", &lastPrime, &ii_Base) != 2)
         FatalError("Line 1 is not a valid newgpen line in input file %s", is_InputTermsFileName.c_str());
      
      ii_N = 0;
                 
      format = FF_NEWPGEN;
      ib_OnlyTwins = true;
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
            if (sscanf(buffer, "%" SCNu64"", &diff) != 1)
               FatalError("Line %s is malformed", buffer);
            
            k += diff;
            break;

         case FF_ABC:
            if (ib_OnlyTwins)
            {
               if (sscanf(buffer, "%" SCNu64"", &k) != 1)
                  FatalError("Line %s is malformed", buffer);
            }
            else
            {
               if (sscanf(buffer, "%" SCNu64" %d", &k, &c) != 2)
                  FatalError("Line %s is malformed", buffer);
            }
            
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
            break;

         default:
            FatalError("Input file %s has unknown format", is_InputTermsFileName.c_str());
      }
            
      if (haveBitMap)
      {
         bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));

         if (ib_OnlyTwins)
         {
            iv_TwinTerms[bit] = true;
            il_TermCount++;
         }
         else
         {
            if (c == -1)
               iv_MinusTerms[bit] = true;
            else
               iv_PlusTerms[bit] = true;
            il_TermCount++;
         }
      }
      else
      {
         if (il_MinK > k) il_MinK = k;
         if (il_MaxK < k) il_MaxK = k;
      }
   }

   fclose(fPtr);
}

void TwinApp::BuildPrimorialTerms(void)
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

void TwinApp::BuildFactorialTerms(void)
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

bool TwinApp::ApplyFactor(uint64_t theFactor,  const char *term)
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
   
   if (c != +1 && c != -1)
      FatalError("Expected c of +1 or -1 in factor but found %d", c);
        
   if (k < il_MinK || k > il_MaxK)
      return false;
   
   uint64_t bit = k - il_MinK;
   
   // No locking is needed because the Workers aren't running yet
   if (ib_OnlyTwins)
   {
      if (iv_TwinTerms[bit])
      {
         iv_TwinTerms[bit] = false;
         il_TermCount--;

         return true;
      }
   }
   else
   {
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
   }
      
   return false;
}

void TwinApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   uint64_t termsCounted = 0;
   
   FILE    *termsFile = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!termsFile)
      FatalError("Unable to open output file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();
      
   if (it_Format == FF_ABCD)
      termsCounted = WriteABCDTermsFile(largestPrime, termsFile);
   
   if (it_Format == FF_ABC)
      termsCounted = WriteABCTermsFile(largestPrime, termsFile);
   
   if (it_Format == FF_NEWPGEN)
      termsCounted = WriteNewPGenTermsFile(largestPrime, termsFile);
   
   fclose(termsFile);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

uint64_t TwinApp::WriteABCDTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0, previousK;
   uint64_t bit, adder = (ib_HalfK ? 2 : 1);

   k = il_MinK;

   bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));

   for (; k<=il_MaxK; k+=adder)
   {      
      if (iv_TwinTerms[bit])
         break;
      
      bit++;
   }

   if (k > il_MaxK)
      return 0;
   
   if (it_TermType == TT_BN)
      fprintf(termsFile, "ABCD $a*%u^%d+1 & $a*%u^%d-1  [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_Base, ii_N, ii_Base, ii_N, k, maxPrime);

   if (it_TermType == TT_PRIMORIAL)
      fprintf(termsFile, "ABCD $a*%u#+1 & $a*%u#-1  [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_N, ii_N, k, maxPrime);

   if (it_TermType == TT_FACTORIAL)
      fprintf(termsFile, "ABCD $a*%u!+1 & $a*%u!-1  [%" SCNu64"] // Sieved to %" SCNu64"\n", ii_N, ii_N, k, maxPrime);
   
   previousK = k;
   kCount = 1;
   k += adder;
   
   bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));
   
   for (; k<=il_MaxK; k+=adder)
   {
      if (iv_TwinTerms[bit])
      {
         fprintf(termsFile, "%" PRIu64"\n", k - previousK);
         previousK = k;
         kCount++;
      }
      
      bit++;
   }

   return kCount;
}

uint64_t TwinApp::WriteABCTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0;
   uint64_t bit, adder = (ib_HalfK ? 2 : 1);

   if (ib_OnlyTwins)
   {
      if (it_TermType == TT_BN)
         fprintf(termsFile, "ABC $a*%u^%d+1 & $a*%u^%d-1 // Sieved to %" SCNu64"\n", ii_Base, ii_N, ii_Base, ii_N, maxPrime);

      if (it_TermType == TT_PRIMORIAL)
         fprintf(termsFile, "ABC $a*%u#+1 & $a*%u#-1 // Sieved to %" SCNu64"\n", ii_N, ii_N, maxPrime);

      if (it_TermType == TT_FACTORIAL)
         fprintf(termsFile, "ABC $a*%u!+1 & $a*%u!-1 // Sieved to %" SCNu64"\n", ii_N, ii_N, maxPrime);
   }
   else
   {
      if (it_TermType == TT_BN)
         fprintf(termsFile, "ABC $a*%u^%d$b // Sieved to %" SCNu64"\n", ii_Base, ii_N, maxPrime);

      if (it_TermType == TT_PRIMORIAL)
         fprintf(termsFile, "ABC $a*%u#$b // Sieved to %" SCNu64"\n", ii_N, maxPrime);

      if (it_TermType == TT_FACTORIAL)
         fprintf(termsFile, "ABC $a*%u!$b // Sieved to %" SCNu64"\n", ii_N, maxPrime);
   }
      
   k = il_MinK;
   bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));

   for ( ; k<=il_MaxK; k+=adder)
   {
      if (ib_OnlyTwins)
      {
         if (iv_TwinTerms[bit])
         {
            fprintf(termsFile, "%" PRIu64"\n", k);
            kCount++;
         }
      }
      else
      {
         if (iv_PlusTerms[bit])
         {
            fprintf(termsFile, "%" PRIu64" +1\n", k);
            kCount++;
         }
         
         if (iv_MinusTerms[bit])
         {
            fprintf(termsFile, "%" PRIu64" -1\n", k);
            kCount++;
         }
      }
      
      bit++;
   }
   
   return kCount;
}

uint64_t TwinApp::WriteNewPGenTermsFile(uint64_t maxPrime, FILE *termsFile)
{
   uint64_t k, kCount = 0;
   uint64_t bit, adder = (ib_HalfK ? 2 : 1);

   fprintf(termsFile, "%" PRIu64":T:0:%u:3\n", maxPrime, ii_Base);
      
   k = il_MinK;
   bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));
   
   for ( ; k<=il_MaxK; k+=adder)
   {
      if (iv_TwinTerms[bit])
      {
         fprintf(termsFile, "%" PRIu64" %u\n", k, ii_N);
         kCount++;
      }
      
      bit++;
   }
   
   return kCount;
}

void  TwinApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{
   if (it_TermType == TT_BN)
      snprintf(extraText, maxTextLength, "%" PRIu64 " < k < %" PRIu64", k*%u^%u", il_MinK, il_MaxK, ii_Base, ii_N);
   
   if (it_TermType == TT_PRIMORIAL)
      snprintf(extraText, maxTextLength, "%" PRIu64" < k < %" PRIu64", k*%u#", il_MinK, il_MaxK, ii_N);
   
   if (it_TermType == TT_FACTORIAL)
      snprintf(extraText, maxTextLength, "%" PRIu64" < k < %" PRIu64", k*%u!", il_MinK, il_MaxK, ii_N);
}

bool  TwinApp::ReportFactor(uint64_t theFactor, uint64_t k, int32_t c)
{
   bool     removedTerm = false;
   char     kStr[50];

   if (theFactor > il_MaxPrimeForValidFactor)
      return false;
   
   // If the first term is valid, then the rest are valid.  In other words 
   // k*x (mod p) = (k+p)*x (mod p) = ... = (k+n*p)*x (mod p)
   VerifyFactor(theFactor, k, c);
   
   ip_FactorAppLock->Lock();

   do
   {
      uint64_t bit = (ib_HalfK ? BIT_HK(k) : BIT_AK(k));

      snprintf(kStr, sizeof(kStr), "%" PRIu64"", k);

      if (ib_OnlyTwins)
      {
         if (iv_TwinTerms[bit])
         {
            iv_TwinTerms[bit] = false;
            removedTerm = true;
            
            LogFactor(theFactor, "%s*%u^%u%+d", kStr, ii_Base, ii_N, c);
            
            il_FactorCount++;
            il_TermCount--;
         }
      }
      else
      {
         if (c == -1 && iv_MinusTerms[bit])
         {
            iv_MinusTerms[bit] = false;
            removedTerm = true;
            
            LogFactor(theFactor, "%s*%u^%u-1", kStr, ii_Base, ii_N);
            
            il_FactorCount++;
            il_TermCount--;
         }

         if (c == +1 && iv_PlusTerms[bit])
         {
            iv_PlusTerms[bit] = false;
            removedTerm = true;
            
            LogFactor(theFactor, "%s*%u^%u+1", kStr, ii_Base, ii_N);
            
            il_FactorCount++;
            il_TermCount--;
         }
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
   
   return removedTerm;
}

// Don't sieve beyond sqrt(maxk*b^n+c).  The worker will not
// remove terms that are prime, so this means that all remaining
// terms will be prime when we reach p = sqrt(maxk*b^n+c).
void  TwinApp::AdjustMaxPrime(void)
{  
   double    toobig = (double) GetMaxPrime();
   double    c = -1.0;
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

void  TwinApp::VerifyFactor(uint64_t theFactor, uint64_t k, int32_t c)
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
   
   if (c == +1)
      resRem = mp.add(resRem, pOne);
   else
      resRem = mp.sub(resRem, pOne);
   
   if (resRem != mp.zero())
   {
      char     term[50];

      if (it_TermType == TT_BN)
         snprintf(term, sizeof(term), "%" PRIu64"*%u^%u%+d", k, ii_Base, ii_N, c);
         
      if (it_TermType == TT_PRIMORIAL)
         snprintf(term, sizeof(term), "%" PRIu64"*%u#%+d", k, ii_N, c);
   
      if (it_TermType == TT_FACTORIAL)
         snprintf(term, sizeof(term), "%" PRIu64"*%u!%+d", k, ii_N, c);

      FatalError("Invalid factor: %" PRIu64" is not a factor of %s", theFactor, term);
   }
}