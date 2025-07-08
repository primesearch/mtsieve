/* HyperCullenWoodallApp.cpp -- (C) Mark Rodenkirch

   HyperCullenWoodallSieve - sieve x^y*y^x+1 or x^y*y^x-1

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <stdint.h>
#include <stdarg.h>
#include <time.h>
#include "../core/inline.h"
#include "../core/MpArith.h"

#include "HyperCullenWoodallApp.h"
#include "HyperCullenWoodallWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "HyperCullenWoodallGpuWorker.h"
#endif

#define APP_NAME        "hcwsieve"
#define APP_VERSION     "1.1"

#define BIT(b, n)       ((((b) - ii_MinB) * GetNCount()) + ((n) - ii_MinN))

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the GPUSieve framework for other applications.
App *get_app(void)
{
   return new HyperCullenWoodallApp();
}

HyperCullenWoodallApp::HyperCullenWoodallApp(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors numbers of the form b^n*n^b+1 and/or b^n*n^b-1");
   SetLogFileName("hcwsieve.log");

   ii_MinB = 0;
   ii_MaxB = 0;
   ii_MinN = 0;
   ii_MaxN = 0;
   ii_CpuWorkSize = 10000;
   ib_IsPlus = false;
   ib_IsMinus = false;
   SetAppMinPrime(3);
   ip_bTerms = NULL;
   ip_nTerms = NULL;
   ib_SplitTerms = false;
   ii_SplitB = 0;
   ii_SplitN = 0;
      
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuSteps = 100000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}
void HyperCullenWoodallApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-b --minb=b           minimum b to search\n");
   printf("-B --maxb=B           maximum b to search\n");
   printf("-n --minn=n           minimum n to search\n");
   printf("-N --maxn=N           maximum n to search\n");
   printf("-s --sign=+/-         sign to sieve for\n");
   printf("-t --splitterms       split terms lower than b and n\n");

#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-S --step=S           max steps iterated per call to GPU (default %d)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

HyperCullenWoodallApp::~HyperCullenWoodallApp(void)
{
    if (ip_bTerms != NULL)
      delete ip_bTerms;
   
    if (ip_nTerms != NULL)
      delete ip_nTerms;
}

void  HyperCullenWoodallApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "b:B:N:n:s:S:M:t";

   AppendLongOpt(longOpts, "MinB",              required_argument, 0, 'b');
   AppendLongOpt(longOpts, "MinN",              required_argument, 0, 'B');
   AppendLongOpt(longOpts, "minn",              required_argument, 0, 'n');
   AppendLongOpt(longOpts, "maxn",              required_argument, 0, 'N');
   AppendLongOpt(longOpts, "sign",              required_argument, 0, 's');
   AppendLongOpt(longOpts, "splitterms",        no_argument, 0, 't');
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   AppendLongOpt(longOpts, "steps",             required_argument, 0, 'S');
   AppendLongOpt(longOpts, "maxfactors",        required_argument, 0, 'M');
#endif
}


parse_t HyperCullenWoodallApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'b':
         status = Parser::Parse(arg, 2, 1000000000, ii_MinB);
         break;

      case 'B':
         status = Parser::Parse(arg, 2, 1000000000, ii_MaxB);
         break;
		 
      case 'n':
         status = Parser::Parse(arg, 3, 1000000000, ii_MinN);
         break;

      case 'N':
         status = Parser::Parse(arg, 3, 1000000000, ii_MaxN);
         break;
         
      case 's':
         char value;
         
         status = Parser::Parse(arg, "+-b", value);
         
         ib_IsPlus = (value == '+' || value == 'b');
         ib_IsMinus = (value == '-' || value == 'b');
         break;

       case 't':
         ib_SplitTerms = true;
         status = P_SUCCESS;
         break;
         
#if defined(USE_OPENCL) || defined(USE_METAL)
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

void HyperCullenWoodallApp::ValidateOptions(void)
{
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "hcw.pfgw";

   if (is_InputTermsFileName.length() > 0)
   {
      if (ib_SplitTerms)
      {
         if (ii_MinB == 0 && ii_MinN == 0)
            FatalError("x and y must be specified if splitting terms");
         
         ii_SplitB = ii_MinB;
         ii_SplitN = ii_MinN;
         ii_MinB = ii_MinN = 0;
      }

      ProcessInputTermsFile(false);

      iv_PlusTerms.resize(GetBCount() * GetNCount());
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
      
      iv_MinusTerms.resize(GetBCount() * GetNCount());
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
         
      il_TermCount = 0;
      
      ProcessInputTermsFile(true);      

      if (ib_SplitTerms)
      {
         WriteOutputTermsFile(il_SplitLargestPrimeTested);
         WriteToConsole(COT_OTHER, "Write %" PRIu64" terms to %s", il_TermCount, is_InputTermsFileName.c_str());
         exit(0);
      }
   }
   else
   {
      if (!ib_IsPlus && !ib_IsMinus)
         FatalError("sign must be specified");

      if (ii_MinB == 0)
         FatalError("min b has not been specified");
	  
      if (ii_MaxB == 0)
         FatalError("max b has not been specified");

      if (ii_MinN == 0)
         FatalError("min n has not been specified");

      if (ii_MaxN == 0)
         FatalError("max n has not been specified");

      if (ii_MinB > ii_MaxB)
         FatalError("min b must be less than or equal to max b");

      if (ii_MinN > ii_MaxN)
         FatalError("min n must be less than or equal to max n");

      if (ii_MinN <= ii_MinB)
         FatalError("min n must be greater than min b");

      if (ii_MaxN <= ii_MaxB)
         FatalError("max n must be greater than max n");         

      il_TermCount = GetBCount() * GetNCount();
      
      if (iv_PlusTerms.max_size() < il_TermCount)
         FatalError("Not enough memory to hold the terms");
      
      iv_PlusTerms.resize(il_TermCount);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
               
      iv_MinusTerms.resize(GetBCount() * GetNCount());
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
      
      SetInitialTerms();
   }
   
   FactorApp::ParentValidateOptions();

   SetMinGpuPrime(100000);
}

Worker *HyperCullenWoodallApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#if defined(USE_OPENCL) || defined(USE_METAL)  
   if (gpuWorker)
      theWorker = new HyperCullenWoodallGpuWorker(id, this);
   else
#endif
      theWorker = new HyperCullenWoodallWorker(id, this);

   return theWorker;
}

void  HyperCullenWoodallApp::SetInitialTerms(void)
{
   uint32_t   b, n;
   uint64_t   bit;
   uint32_t   xgtxCount = 0;
   uint32_t   evenCount = 0;
   uint32_t   algebraicCount = 0;
   
   // Reset this
   il_TermCount = 0;

   for (b=ii_MinB; b<=ii_MaxB; b++)
   {
      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         bool includeMinus = true;

         if (n <= b)
         {
            if (ib_IsPlus) xgtxCount++;
            if (ib_IsMinus) xgtxCount++;
            continue;
         }
         
         // If b and n are both odd, then all terms are even
         if (b & 1 && n & 1)
         {
            if (ib_IsPlus) evenCount++;
            if (ib_IsMinus) evenCount++;
            continue;
         }
         
         // If b and n are both even, then -1 terms have an algebraic factor
         if (!(b & 1) && !(n & 1))
         {
            if (ib_IsMinus)
            {
               algebraicCount++;
               includeMinus = false;
            }
         }

         bit = BIT(b, n);
   
         iv_PlusTerms[bit] = true;            
         il_TermCount++;
         
         if (includeMinus)
         {
            iv_MinusTerms[bit] = true;            
            il_TermCount++;
         }
      }
   }

   WriteToConsole(COT_OTHER, "Quick elimination of terms info (in order of check):");
   WriteToConsole(COT_OTHER, "    %u because n <= b", xgtxCount);
   WriteToConsole(COT_OTHER, "    %u because the term is even", evenCount);
   WriteToConsole(COT_OTHER, "    %u because the term has an algebraic factor", algebraicCount);
}

void HyperCullenWoodallApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   FILE    *stPtr = NULL;
   char     buffer[1000];
   uint32_t b, n;
   int32_t  sign;
   uint64_t minPrime;
   uint32_t splitCount = 0;
   
   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("File %s is empty", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC $a^$b*$b^$a$c // Sieved to %" SCNu64"", &minPrime) != 1)
      FatalError("First line of the input file is malformed");

   if (!haveBitMap)
   {
      ib_IsPlus = false;
      ib_IsMinus = false;
   }

   // Reset this
   il_TermCount = 0;

   SetMinPrime(minPrime);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (sscanf(buffer, "%u %u %d", &b, &n, &sign) != 3)
         FatalError("Line %s is malformed", buffer);

      if (sign != +1 && sign != -1)
         FatalError("Only +1 or -1 is supported");

      if (ib_SplitTerms && b <= ii_SplitB && n <= ii_SplitN)
      {
         if (!haveBitMap)
         {
            if (stPtr == NULL)
            {
               stPtr = fopen("split_hcw.pfgw", "w");
               fprintf(stPtr, "ABC $a^$b*$b^$a$c // Sieved to %" PRIu64"\n", minPrime);         
            }
            
            fprintf(stPtr, "%s\n", buffer);
            splitCount++;
         }
         continue;
      }

      if (haveBitMap)
      {
         if (sign == +1)
            iv_PlusTerms[BIT(b, n)] = true;
         else
            iv_MinusTerms[BIT(b, n)] = true;
      }
      else
      {
         if (sign == +1)
            ib_IsPlus = true;
         
         if (sign == -1)
            ib_IsMinus = true;
         
         if (ii_MinB == 0 || b < ii_MinB)
            ii_MinB = b;

         if (b > ii_MaxB)
            ii_MaxB = b;

         if (ii_MinN == 0 || n < ii_MinN)
            ii_MinN = n;

         if (n > ii_MaxN)
            ii_MaxN = n;
      }
         
      il_TermCount++;
   }

   fclose(fPtr);

   if (ib_SplitTerms)
   {
      if (splitCount > 0)
      {
         fclose(stPtr);
         WriteToConsole(COT_OTHER, "Split off %u terms to split_hcw.pfgw", splitCount);
      }

      il_SplitLargestPrimeTested = minPrime;
   }
}

bool HyperCullenWoodallApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint32_t b1, b2, n1, n2;
   int32_t  sign;
   
   if (sscanf(term, "%u^%u*%u^%u%d", &b1, &n1, &n2, &b2, &sign) != 5)
      FatalError("Could not parse term %s\n", term);

   if (b1 != b2)
      FatalError("x values for term %s do not match (%u != %u)\n", term, b1, b2);
   
   if (n1 != n2)
      FatalError("y values for term %s do not match (%u != %u)\n", term, n1, n2);

   if (sign != +1 && sign != -1)
      FatalError("sign term %s must be +1 or -1\n", term);
   
   if (b1 < ii_MinB || b1 > ii_MaxB)
      return false;
   
   if (n1 < ii_MinN || n1 > ii_MaxN)
      return false;

   if (!ib_IsMinus && sign != +1)
      FatalError("+/- for term %s is not the expected value\n", term);
   
   if (!ib_IsPlus && sign != -1)
      FatalError("+/- for term %s is not the expected value\n", term);

   VerifyFactor(theFactor, b1, n1, sign);

   uint64_t bit = BIT(b1, n1);
   bool haveTerm = false;
   
   // No locking is needed because the Workers aren't running yet
   if (sign > 0)
      haveTerm = iv_PlusTerms[bit];
   else
      haveTerm = iv_MinusTerms[bit];

   if (haveTerm)
   {
      if (sign > 0)
         iv_PlusTerms[bit] = false;
      else
         iv_MinusTerms[bit] = false;
      
      il_TermCount--;
      return true;
   }

   return false;
}

void HyperCullenWoodallApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *fPtr;
   uint32_t b, n;
   uint64_t bit, termsCounted = 0;
   
   ip_FactorAppLock->Lock();

   fPtr = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!fPtr)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   fprintf(fPtr, "ABC $a^$b*$b^$a$c // Sieved to %" PRIu64"\n", largestPrime);

   for (b=ii_MinB; b<=ii_MaxB; b++)
   {
      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         bit = BIT(b, n);
         
         if (iv_PlusTerms[bit])
         {
            fprintf(fPtr, "%u %u +1\n", b, n);
            termsCounted++;
         }
         
         if (iv_MinusTerms[bit])
         {
            fprintf(fPtr, "%u %u -1\n", b, n);
            termsCounted++;
         }
      }
   }

   fclose(fPtr);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);

   ip_FactorAppLock->Release();
}

void  HyperCullenWoodallApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{   
   if (!ib_IsMinus)
      snprintf(extraText, maxTextLength, "%d <= b <= %d, %d <= n <= %d (+1 only)",ii_MinB, ii_MaxB, ii_MinN, ii_MaxN);
   else if (!ib_IsPlus)
      snprintf(extraText, maxTextLength, "%d <= b <= %d, %d <= n <= %d (-1 only)",ii_MinB, ii_MaxB, ii_MinN, ii_MaxN);
   else
      snprintf(extraText, maxTextLength, "%d <= b <= %d, %d <= n <= %d (both +1 and -1)",ii_MinB, ii_MaxB, ii_MinN, ii_MaxN);
}

bool HyperCullenWoodallApp::ReportFactor(uint64_t theFactor, uint32_t b, uint32_t n, int32_t sign)
{
   uint64_t bit;
   bool     removedTerm = false;
   
   if (b < ii_MinB || b > ii_MaxB)
      return false;
   
   if (n < ii_MinN || n > ii_MaxN)
      return false;
   
   if (sign != -1 && sign != +1)
      FatalError("expecting sign of +1 or -1");
   
   if (sign == +1 && !ib_IsPlus)
      return false;
   
   if (sign == -1 && !ib_IsMinus)
      return false;
   
   VerifyFactor(theFactor, b, n, sign);
   
   ip_FactorAppLock->Lock();

   bit = BIT(b, n);
   
   bool haveTerm = false;

   if (ib_IsPlus && sign == +1)
      haveTerm = iv_PlusTerms[bit];
   
   if (ib_IsMinus && sign == -1)
      haveTerm = iv_MinusTerms[bit];
   
   if (haveTerm)
   {
      if (sign > 0)
         iv_PlusTerms[bit] = false;
      else
         iv_MinusTerms[bit] = false;

      il_TermCount--;
      il_FactorCount++;
      removedTerm = true;
      LogFactor(theFactor, "%u^%u*%u^%u%+d", b, n, n, b, sign);
   }
   
   ip_FactorAppLock->Release();

   return removedTerm;
}

void  HyperCullenWoodallApp::VerifyFactor(uint64_t theFactor, uint32_t b, uint32_t n, int32_t sign)
{
   bool     isValid = true;
   
   MpArith  mp(theFactor);
   MpRes    resBN = mp.pow(mp.nToRes(b), n);
   MpRes    resNB = mp.pow(mp.nToRes(n), b);
   
   resBN = mp.mul(resBN, resNB);

   uint64_t res = mp.resToN(resBN);
   
   if (sign > 0 && res != theFactor - 1)
      isValid = false;
   
   if (sign < 0 && res != 1)
      isValid = false;
   
   if (isValid)
      return;
   
   FatalError("%" PRIu64" does not divide %u^%u*%u^%u%+d", theFactor, b, n, n, b, sign);
}

void  HyperCullenWoodallApp::BuildTerms(void)
{
   uint32_t bit;

   // This could be done using a lot less memory, but the time and effort to do that
   // is significant.  That will be left as a future task for anyone intested in doing it.
   ip_FactorAppLock->Lock();

   if (ip_bTerms == NULL)
      ip_bTerms = (base_t *) xmalloc(ii_MaxB - ii_MinB + 2, sizeof(base_t) , "bTerms");

   base_t *bTerm = ip_bTerms;

   for (uint32_t b=ii_MinB; b<=ii_MaxB; b++)
   {      
      bTerm->base = b;
      bTerm->powerCount = 0;
      bTerm->indexedByPower = false;

      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         bit = BIT(b, n);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            bTerm->powerCount++;
      }

      bTerm->powers.resize(bTerm->powerCount);
      bTerm->residues.resize(bTerm->powerCount);
      
      bTerm->powerCount = 0;

      for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
      {
         bit = BIT(b, n);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            bTerm->powers[bTerm->powerCount] = n;
            bTerm->powerCount++;
         }
      }

      bTerm++;
   }

   if (ip_nTerms == NULL)
      ip_nTerms = (base_t *) xmalloc(ii_MaxN - ii_MinN + 2, sizeof(base_t) , "nTerms");
   
   base_t *nTerm = ip_nTerms;
   
   for (uint32_t n=ii_MinN; n<=ii_MaxN; n++)
   {
      nTerm->base = n;

      for (uint32_t b=ii_MinB; b<=ii_MaxB; b++)
      {
         bit = BIT(b, n);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            nTerm->powerCount++;
      }

      nTerm->powers.resize(nTerm->powerCount);
      nTerm->residues.resize(ii_MaxB - ii_MinB + 1);
      nTerm->indexedByPower = true;

      nTerm->powerCount = 0;

      for (uint32_t b=ii_MinB; b<=ii_MaxB; b++)
      {
         bit = BIT(b, n);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            nTerm->powers[nTerm->powerCount] = b;
            nTerm->powerCount++;
         }
      }
      
      nTerm++;
   }

   ip_FactorAppLock->Release();
}

#if defined(USE_OPENCL) || defined(USE_METAL)
void     HyperCullenWoodallApp::FreeGpuGroupsOfTerms(gpugroup_t *gPtr)
{
   gpugroup_t *cgPtr = NULL;

   while (gPtr != NULL)
   {
      cgPtr = gPtr;
      gPtr = (gpugroup_t *) gPtr->nextGroup;
      
      xfree(cgPtr->terms);
      xfree(cgPtr);
   }
}
   
// Build terms in groups.  Each group will have:
//   b1 b1n1 b1n2 ... b1yn 0
//   b2 b2n1 b2n2 ... b2yn 0
//    ... ... ... ... ...
//   xm xmn1 xmn2 ... xmyn 0
//   0
// Each group of terms for x will end with a 0 indicating that the next term
// after the 0 is for the next x.  Two 0s in a row indicates end of the group.
gpugroup_t  *HyperCullenWoodallApp::GetGpuGroupsOfTerms(void)
{
   uint32_t    bit, b, n;
   uint32_t    termsForB = 0;
   gpugroup_t *fgPtr = NULL;
   gpugroup_t *cgPtr = NULL;
   uint32_t    index;
   uint32_t    distinctBN = 0;
   uint32_t   *terms;
   int         groups = 1;
   
   ip_FactorAppLock->Lock();

   index = 0;

   b = ii_MinB;
   
   fgPtr = (gpugroup_t *) xmalloc(1, sizeof(gpugroup_t), "terms");
   fgPtr->terms = (uint32_t *) xmalloc(ii_MaxGpuSteps, sizeof(uint32_t), "steps");
   fgPtr->nextGroup = NULL;

   cgPtr = fgPtr;

   terms = fgPtr->terms;
   do
   {
      termsForB = 0;
      
      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         bit = BIT(b, n);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            termsForB++;
      }

      if (termsForB == 0)
         continue;

      // If not enough space for all y for this x, then put this x into the next group.
      if ((termsForB + index) >= (ii_MaxGpuSteps - 2))
      {
         if (index == 0)
            FatalError("-S must be at least %u", termsForB+2);

         // Signal that there are no more x for this group
         terms[index] = 0;
         
         cgPtr->nextGroup = (gpugroup_t *) xmalloc(1, sizeof(gpugroup_t), "terms");
         
         cgPtr = (gpugroup_t *) cgPtr->nextGroup;
         cgPtr->terms = (uint32_t *) xmalloc(ii_MaxGpuSteps, sizeof(uint32_t), "steps");
         index = 0;
         
         terms = cgPtr->terms;
         groups++;
      }
      
      terms[index] = b;
      index++;

      for (n=ii_MinN; n<=ii_MaxN; n++)
      {
         bit = BIT(b, n);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            terms[index] = n;
            index++;
            distinctBN++;
         }
      }

      // Signal that there are no more n for this b
      terms[index] = 0;
      index++;
      
      b++;
   } while (b <= ii_MaxB);

   ip_FactorAppLock->Release();

   WriteToConsole(COT_OTHER, "%u GPU group%s needed to hold %u distinct b/n for %llu terms", groups, (groups > 1 ? "s are" : " is"), distinctBN, il_TermCount);
   return fgPtr;
}
#endif
