/* CullenWoodallApp.cpp -- (C) Mark Rodenkirch, September 2012

   This program finds factors of numbers in the form n*b^n+1 and n*b^n-1

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <stdint.h>
#include <stdarg.h>
#include <time.h>
#include "../x86_asm/fpu-asm-x86.h"

#include "CullenWoodallApp.h"
#include "CullenWoodallWorker.h"

#ifdef HAVE_GPU_WORKERS
#include "CullenWoodallGpuWorker.h"
#define APP_NAME        "gcwsievecl"
#else
#define APP_NAME        "gcwsieve"
#endif

#define APP_VERSION     "1.4"

#define BIT(n)        ((n) - ii_MinN)

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the GPUSieve framework for other applications.
App *get_app(void)
{
   return new CullenWoodallApp();
}

CullenWoodallApp::CullenWoodallApp(void) : AlgebraicFactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors numbers of the form n*b^n+1 and n*b^n-1");
   SetLogFileName("gcwsieve.log");

   ii_Base = 0;
   ii_MinN = 0;
   ii_MaxN = 0;
   ib_Cullen = false;
   ib_Woodall = false;
   it_Format = FF_ABC;

   SetAppMinPrime(3);
   
#ifdef HAVE_GPU_WORKERS
   ii_MaxGpuSteps = 100000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 10;
#endif
}

void CullenWoodallApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-b --base=b           Base to search\n");
   printf("-n --min_n=n          Minimum n to search\n");
   printf("-N --max_n=N          Maximum N to search\n");
   printf("-s --sign=+/-/b       Sign to sieve for (+ = Cullen, - = Woodall)\n");
   printf("-f --format=f         Format of output file (A=ABC (default), L=LLR\n");
#ifdef HAVE_GPU_WORKERS
   printf("-S --step=S           max steps iterated per call to GPU (default %d)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  CullenWoodallApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "b:n:N:s:";

   AppendLongOpt(longOpts, "base",           required_argument, 0, 'b');
   AppendLongOpt(longOpts, "min_n",          required_argument, 0, 'n');
   AppendLongOpt(longOpts, "max_n",          required_argument, 0, 'N');
   AppendLongOpt(longOpts, "sign",           required_argument, 0, 's');
   AppendLongOpt(longOpts, "format",         required_argument, 0, 'f');
   
#ifdef HAVE_GPU_WORKERS
   shortOpts += "S:M:";
   
   AppendLongOpt(longOpts, "maxsteps",       required_argument, 0, 'S');
   AppendLongOpt(longOpts, "maxfactors",     required_argument, 0, 'M');
#endif
}


parse_t CullenWoodallApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;
   char value;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'b':
         status = Parser::Parse(arg, 2, 1000000000, ii_Base);
         break;

      case 'n':
         status = Parser::Parse(arg, 2, 1000000000, ii_MinN);
         break;

      case 'N':
         status = Parser::Parse(arg, 2, 1000000000, ii_MaxN);
         break;
		 
      case 'f':
         status = Parser::Parse(arg, "AL", value);
         
         it_Format = FF_UNKNOWN;
   
         if (value == 'A')
            it_Format = FF_ABC;
         if (value == 'L')
            it_Format = FF_LLR;
         break;
         
      case 's':
         status = Parser::Parse(arg, "+-b", value);
         if (value == '-')
            ib_Woodall = true;
         if (value == '+')
            ib_Cullen = true;
         if (value == 'b')
            ib_Woodall = ib_Cullen = true;
         break;
         
#ifdef HAVE_GPU_WORKERS
      case 'S':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxGpuSteps);
         break;

      case 'M':
         status = Parser::Parse(arg, 10, 1000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void CullenWoodallApp::ValidateOptions(void)
{
   if (it_Format == FF_UNKNOWN)
      FatalError("the specified file format in not valid, use A (ABC) or L (LLR)");
      
   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      il_TermCount = ii_MaxN - ii_MinN + 1;
      
      iv_CullenTerms.resize(il_TermCount);
      std::fill(iv_CullenTerms.begin(), iv_CullenTerms.end(), false);

      iv_WoodallTerms.resize(il_TermCount);
      std::fill(iv_WoodallTerms.begin(), iv_WoodallTerms.end(), false);

      il_TermCount = 0;
      
      ProcessInputTermsFile(true);
   }
   else
   {
      if (ii_Base == 0)
         FatalError("base has not been specified");

      if (ii_MinN == 0)
         FatalError("min n has not been specified");

      if (ii_MaxN == 0)
         FatalError("max n has not been specified");

      if (ii_MinN > ii_MaxN)
         FatalError("min n cannot be greater than max n");

      if (!ib_Cullen && !ib_Woodall)
         FatalError("Choose Cullen and/or Woodall form");

      il_TermCount = ii_MaxN - ii_MinN + 1;
      
      iv_CullenTerms.resize(il_TermCount);
      std::fill(iv_CullenTerms.begin(), iv_CullenTerms.end(), false);
      
      iv_WoodallTerms.resize(il_TermCount);
      std::fill(iv_WoodallTerms.begin(), iv_WoodallTerms.end(), false);
      
      SetInitialTerms();
      
      EliminateGfnAndMersenneTerms();
   }
      
   if (is_OutputTermsFileName.length() == 0)
   {
      char fileName[20];
      
      sprintf(fileName, "gcw_b%u.pfgw", ii_Base);
      is_OutputTermsFileName = fileName;
   }
   
   FactorApp::ParentValidateOptions();
   
   // This is only done when starting a new sieve, but has to be done
   // after the call to ParentValidateOptions() as that method will 
   // open the factor file that algebraic factors will be written to.
   if (is_InputTermsFileName.length() == 0)
      EliminateTermsWithAlgebraicFactors();

   uint32_t pForSingleWorker = ((ii_Base < ii_MaxN) ? (ii_MaxN + 1) : (ii_Base + 1));
   
   // The Worker will trigger a rebuild of terms when it reaches this prime
   // At this prime multiple threads can be used.
   SetMaxPrimeForSingleWorker(pForSingleWorker);
   
   // The GPU code for this sieve will not support primes lower than this.
   SetMinGpuPrime(pForSingleWorker);
}

Worker *CullenWoodallApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#ifdef HAVE_GPU_WORKERS  
   if (gpuWorker)
      theWorker = new CullenWoodallGpuWorker(id, this);
   else
#endif
      theWorker = new CullenWoodallWorker(id, this);
   
   return theWorker;
}

void CullenWoodallApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t n, b;
   int32_t  c;
   uint64_t p;
   format_t format = FF_UNKNOWN;

   ii_Base = 0;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("File %s is empty", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC $a*%d^$a$b // Sieved to %" SCNu64"", &ii_Base, &p) == 2)
      format = FF_ABC;
      
   if (sscanf(buffer, "ABC $a*$b^$a$c // Sieved to %" SCNu64"", &p) == 1)
      format = FF_LLR;

   if (format == FF_UNKNOWN)
      FatalError("First line of the input file is malformed");

   // Reset this
   il_TermCount = 0;

   SetMinPrime(p);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (format == FF_ABC)
      {
         if (sscanf(buffer, "%u %d", &n, &c) != 2)
            FatalError("Line %s is malformed", buffer);
      }
      
      if (format == FF_LLR)
      {
         if (sscanf(buffer, "%u %u %d", &n, &b, &c) != 3)
            FatalError("Line %s is malformed", buffer);
         
         if (il_TermCount == 0)
            ii_Base = b;
         
         if (ii_Base != b)
            FatalError("Multiple bases specified in input file", buffer);
      }
      

      if (haveBitMap)
      {
         if (c == +1)
            iv_CullenTerms[BIT(n)] = true;
         
         if (c == -1)
            iv_WoodallTerms[BIT(n)] = true;
            
         il_TermCount++;
      }
      else
      {
         if (ii_MinN == 0 || n < ii_MinN)
            ii_MinN = n;

         if (n > ii_MaxN)
            ii_MaxN = n;

         if (c == +1)
            ib_Cullen = true;

         if (c == -1)
            ib_Woodall = true;
      }
   }

   fclose(fPtr);
}

bool CullenWoodallApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t n1, b, n2;
   int32_t  c;
   
   if (sscanf(term, "%u*%u^%u%d", &n1, &b, &n2, &c) != 4)
      FatalError("Could not parse term %s\n", term);

   if (b != ii_Base)
      FatalError("base is correct for term %s (%u != %u)\n", term, ii_Base, b);
   
   if (c != -1 && c != +1)
      FatalError("c is neither +1 nor -1 for term %s\n", c, term);
   
   if (n1 != n2)
      FatalError("n values for term %s do not match (%u != %u)\n", term, n1, n2);
   
   if (n1 < ii_MinN || n1 > ii_MaxN)
      return false;
   
   if (!VerifyFactor(false, thePrime, n1, c))
      return false;
      
   uint64_t bit = BIT(n1);
   
   // No locking is needed because the Workers aren't running yet
   if (c == +1 && iv_CullenTerms[bit])
   {
      iv_CullenTerms[bit] = false;
      il_TermCount--;
      return true;
   }
   
   if (c == -1 && iv_WoodallTerms[bit])
   {
      iv_WoodallTerms[bit] = false;
      il_TermCount--;
      return true;
   }

   return false;
}

void CullenWoodallApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *fPtr = fopen(is_OutputTermsFileName.c_str(), "w");
   uint32_t n, bit;
   uint64_t terms = 0;

   if (!fPtr)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());

   ip_FactorAppLock->Lock();

   if (it_Format == FF_ABC)
      fprintf(fPtr, "ABC $a*%u^$a$b // Sieved to %" PRIu64"\n", ii_Base, largestPrime);
   else
      fprintf(fPtr, "ABC $a*$b^$a$c // Sieved to %" PRIu64"\n", largestPrime);

   for (n=ii_MinN; n<=ii_MaxN; n++)
   {
      bit = BIT(n);
      
      if (iv_CullenTerms[bit])
      {
         terms++;
         
         if (it_Format == FF_ABC)
            fprintf(fPtr, "%u +1\n", n);
         else
            fprintf(fPtr, "%u %u +1\n", n, ii_Base);
      }
      
      if (iv_WoodallTerms[bit])
      {
         terms++;
         
         if (it_Format == FF_ABC)
            fprintf(fPtr, "%u -1\n", n);
         else
            fprintf(fPtr, "%u %u -1\n", n, ii_Base);
      }
   }

   fclose(fPtr);

   if (terms != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%u) != expected terms (%u)", terms, il_TermCount);
   
   ip_FactorAppLock->Release();
}

void CullenWoodallApp::GetExtraTextForSieveStartedMessage(char *extraTtext)
{
   if (ib_Cullen && ib_Woodall)
      sprintf(extraTtext, "%u <= n <= %u, n*%u^n+/-1", ii_MinN, ii_MaxN, ii_Base);
   else if (ib_Cullen)
      sprintf(extraTtext, "%u <= n <= %u, n*%u^n+1", ii_MinN, ii_MaxN, ii_Base);
   else
      sprintf(extraTtext, "%u <= n <= %u, n*%u^n-1", ii_MinN, ii_MaxN, ii_Base);
}

void  CullenWoodallApp::NotifyAppToRebuild(void)
{
}

void  CullenWoodallApp::SetInitialTerms(void)
{
   uint32_t   n;
   uint32_t   evenCount = 0;

   // Reset this
   il_TermCount = 0;

   for (n=ii_MinN; n<=ii_MaxN; n++)
   {
      // If n and b are both odd, then both n*b^n+1 and n*b^n-1 are even
      // so we don't need to add them
      if (n & 1 && ii_Base & 1)
      {
         if (ib_Cullen) evenCount++;
         if (ib_Woodall) evenCount++;
         continue;
      }

      if (ib_Cullen)
      {
         il_TermCount++;
         iv_CullenTerms[BIT(n)] = true;
      }

      if (ib_Woodall)
      {
         il_TermCount++;
         iv_WoodallTerms[BIT(n)] = true;
      }
   }

   if (evenCount > 0)
      WriteToConsole(COT_OTHER, "%d terms removed because the term is even", evenCount);
}

void  CullenWoodallApp::EliminateGfnAndMersenneTerms(void)
{
   uint32_t n, bit;
   uint32_t removedCount = 0;
   
   for (n=ii_MinN; n<=ii_MaxN; n++)
   {
      bit = BIT(n);
      
      if (iv_CullenTerms[bit] && IsGfnOrMersenneForm(n, ii_Base, +1))
      {
         il_TermCount--;
         iv_CullenTerms[bit] = false;
         removedCount++;
      }
      
      if (iv_WoodallTerms[bit] && IsGfnOrMersenneForm(n, ii_Base, -1))
      {
         il_TermCount--;
         iv_WoodallTerms[bit] = false;
         removedCount++;
      }
   }
   
   if (removedCount > 0)
      WriteToConsole(COT_OTHER, "%u terms removed as they were of GFN or Mersenne form", removedCount);
}

// Generalized Woodalls and Cullens might have some easy to find large factors.
// This function will find factors based upon the following:
//    Given a*b^a-1, if a=c^x*b^y and (a+y)%x = 0 for integers c, x, and y,
//    then a*b^a-1 has a factor of c*b^((a+y)/x)-1
//    Given a*b^a+1, if a=c^x*b^y and (a+y)%x = 0 and x is odd for integers
//    c, x, and y, a*b^a+1 has a factor of c*b^((a+y)/x)+1
void  CullenWoodallApp::EliminateTermsWithAlgebraicFactors(void)
{
   uint32_t  b, n, nbase, bpow, npow;
   uint64_t  bexp, nexp;
   uint32_t  removedCount = 0;

   b = ii_Base;

   // Start with nexp=n^0, going to n^1, n^2, etc.
   for (bpow=0, bexp=1; bexp*bexp<=ii_MaxN; bpow++, bexp*=b)
   {
      // Start with nbase=2, going to 3, 4, etc.
      for (nbase=2; nbase*nbase*bexp<=ii_MaxN; nbase++)
      {
         // Start with nexp=nbase^2, going to nbase^3, etc.
         for (npow=2, nexp=nbase*nbase; nexp*bexp<=ii_MaxN; npow++, nexp*=nbase)
         {
            n = nexp * bexp;
            
            if ((n+bpow)%npow == 0)
            {
               if (npow%2 == 1)
                  removedCount += CheckAlgebraicFactor(n, +1, "%u*%u^%u+1", nbase, b, (n+bpow)/npow);
               
               removedCount += CheckAlgebraicFactor(n, -1, "%u*%u^%u-1", nbase, b, (n+bpow)/npow);
            }
         }
      }
   }
   
   if (removedCount > 0)
      WriteToConsole(COT_OTHER, "%u terms removed due to algebraic factorizations", removedCount);
}

bool  CullenWoodallApp::CheckAlgebraicFactor(uint32_t n, int32_t c, const char *fmt, ...)
{
   va_list args;
   char    factor[200];
   bool    removedTerm = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   if (c != +1 && c != -1)
      return false;
   
   va_start(args,fmt);
   vsprintf(factor, fmt, args);
   va_end(args);
      
   uint32_t bit = BIT(n);
   
   if (c == +1 && iv_CullenTerms[bit])
   {
      iv_CullenTerms[bit] = false;
      removedTerm = true;
      
      LogFactor(factor, "%u*%u^%u+1", n, ii_Base, n);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   if (c == -1 && iv_WoodallTerms[bit])
   {
      iv_WoodallTerms[bit] = false;
      removedTerm = true;
      
      LogFactor(factor, "%u*%u^%u-1", n, ii_Base, n);
      
      il_FactorCount++;
      il_TermCount--;
   }
   
   return removedTerm;
}

// Build a one dimensional array containing the values for n that don't have
// a factor.  The 0 at the end indicates that they are no more terms.  This
// array will have n in descending order.
uint32_t  CullenWoodallApp::GetTerms(uint32_t *terms, uint32_t maxTermsInGroup, uint32_t groupSize)
{
   uint32_t index = 0, bit;
   uint32_t groupCount = 0;
   uint32_t termsInGroup = 0, n;
   
   ip_FactorAppLock->Lock();

   for (n=ii_MaxN; n>=ii_MinN; n--)
   {
      // If the group is full, then start another group
      if (termsInGroup > maxTermsInGroup)
      {
         // Indicate that there are no more terms for this group
         terms[index] = 0;
         groupCount++;
         
         index = groupSize * groupCount;
         termsInGroup = 0;
      }

      bit = BIT(n);
      
      if (iv_CullenTerms[bit] || iv_WoodallTerms[bit])
      {
         terms[index] = n;
         index++;
         termsInGroup++;
      }
   }
   
   // Indicate that there are no more terms for the last group
   terms[index] = 0;

   ip_FactorAppLock->Release();

   return groupCount + 1;
}

bool CullenWoodallApp::ReportFactor(uint64_t p, uint32_t n, int32_t c)
{
   uint64_t bit;
   bool     removedTerm = false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
      
   ip_FactorAppLock->Lock();

   bit = BIT(n);
   
   if (ib_Cullen && c == +1 && iv_CullenTerms[bit])
   {
      iv_CullenTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      removedTerm = true;
      LogFactor(p, "%u*%u^%u+1", n, ii_Base, n);

      VerifyFactor(true, p, n, c);
   }
   
   if (ib_Woodall && c == -1 && iv_WoodallTerms[bit])
   {
      iv_WoodallTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      removedTerm = true;
      LogFactor(p, "%u*%u^%u-1", n, ii_Base, n);

      VerifyFactor(true, p, n, c);
   }
   
   ip_FactorAppLock->Release();
   
   return removedTerm;
}

bool  CullenWoodallApp::VerifyFactor(bool badFactorIsFatal, uint64_t p, uint32_t n, int32_t c)
{
   uint64_t rem;
   bool     isValid = true;
      
   fpu_push_1divp(p);
   
   rem = fpu_powmod(ii_Base, n, p);
   rem = fpu_mulmod(rem, n, p);
   
   fpu_pop();
   
   if (c == -1)
      isValid = (rem == +1);
      
   if (c == +1)
      isValid = (rem == p-1);
      
   if (isValid)
      return isValid;
   
   char buffer[100];
   
   sprintf(buffer, "%" PRIu64" is not a factor of %u*%u^%u%+d", p, n, ii_Base, n, c);
   
   if (badFactorIsFatal)
      FatalError(buffer);
   else
      WriteToConsole(COT_OTHER, buffer);

   return isValid;
}
