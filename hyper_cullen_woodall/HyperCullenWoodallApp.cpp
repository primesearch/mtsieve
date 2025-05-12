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
#define APP_VERSION     "1.0"

#define BIT(x, y)       ((((x) - ii_MinX) * GetYCount()) + ((y) - ii_MinY))

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the GPUSieve framework for other applications.
App *get_app(void)
{
   return new HyperCullenWoodallApp();
}

HyperCullenWoodallApp::HyperCullenWoodallApp(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors numbers of the form x^y*y^x+1 and/or x^y*y^x-1");
   SetLogFileName("hcwsieve.log");

   ii_MinX = 0;
   ii_MaxX = 0;
   ii_MinY = 0;
   ii_MaxY = 0;
   ii_CpuWorkSize = 10000;
   ib_IsPlus = false;
   ib_IsMinus = false;
   SetAppMinPrime(3);
   ip_xTerms = NULL;
   ip_yTerms = NULL;
      
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_MaxGpuSteps = 100000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 100;
#endif
}
void HyperCullenWoodallApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-x --minx=x           minimum x to search\n");
   printf("-X --maxx=X           maximum x to search\n");
   printf("-y --miny=y           minimum y to search\n");
   printf("-Y --maxy=Y           maximum y to search\n");
   printf("-s --sign=+/-         sign to sieve for\n");
#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-S --step=S           max steps iterated per call to GPU (default %d)\n", ii_MaxGpuSteps);
   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

HyperCullenWoodallApp::~HyperCullenWoodallApp(void)
{
    if (ip_xTerms != NULL)
      delete ip_xTerms;
   
    if (ip_yTerms != NULL)
      delete ip_yTerms;
}

void  HyperCullenWoodallApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "x:X:y:Y:s:S:M:";

   AppendLongOpt(longOpts, "minx",              required_argument, 0, 'x');
   AppendLongOpt(longOpts, "maxx",              required_argument, 0, 'X');
   AppendLongOpt(longOpts, "miny",              required_argument, 0, 'y');
   AppendLongOpt(longOpts, "sign",              required_argument, 0, 's');
   AppendLongOpt(longOpts, "sign",              required_argument, 0, 's');
   
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
      case 'x':
         status = Parser::Parse(arg, 2, 1000000000, ii_MinX);
         break;

      case 'X':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxX);
         break;
		 
      case 'y':
         status = Parser::Parse(arg, 2, 1000000000, ii_MinY);
         break;

      case 'Y':
         status = Parser::Parse(arg, 1, 1000000000, ii_MaxY);
         break;
         
      case 's':
         char value;
         
         status = Parser::Parse(arg, "+-b", value);
         
         ib_IsPlus = (value == '+' || value == 'b');
         ib_IsMinus = (value == '-' || value == 'b');
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
      ProcessInputTermsFile(false);

      iv_PlusTerms.resize(GetXCount() * GetYCount());
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
      
      iv_MinusTerms.resize(GetXCount() * GetYCount());
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
         
      il_TermCount = 0;
      
      ProcessInputTermsFile(true);      
   }
   else
   {
      if (!ib_IsPlus && !ib_IsMinus)
         FatalError("sign must be specified");

      if (ii_MinX == 0)
         FatalError("min x has not been specified");
	  
      if (ii_MaxX == 0)
         FatalError("max x has not been specified");

      if (ii_MinY == 0)
         FatalError("min y has not been specified");

      if (ii_MaxY == 0)
         FatalError("max y has not been specified");

      if (ii_MinX > ii_MaxX)
         FatalError("min x must be less than or equal to max x");

      if (ii_MinY > ii_MaxY)
         FatalError("min y must be less than or equal to max y");

      il_TermCount = GetXCount() * GetYCount();
      
      if (iv_PlusTerms.max_size() < il_TermCount)
         FatalError("Not enough memory to hold the terms");
      
      iv_PlusTerms.resize(il_TermCount);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
               
      iv_MinusTerms.resize(GetXCount() * GetYCount());
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
   uint32_t   x, y;
   uint64_t   bit;
   uint32_t   evenCount = 0;
   uint32_t   algebraicCount = 0;
   
   // Reset this
   il_TermCount = 0;

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         boolean includeMinus = true;

         // If x and y are both odd, then both x^y*y^x+1 and x^y*y^x-1 are even
         // so we don't need to add them
         if (x & 1 && y & 1)
         {
            if (ib_IsPlus) evenCount++;
            if (ib_IsMinus) evenCount++;
            continue;
         }
         
         // If x and y are both even, then x^y*y^x-1 has an algebraic factor
         // so we don't need to add them
         if (!(x & 1) && !(y & 1))
         {
            if (ib_IsMinus)
            {
               algebraicCount++;
               includeMinus = false;
            }
         }

         bit = BIT(x, y);
   
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
   WriteToConsole(COT_OTHER, "    %u because the term is even", evenCount);
   WriteToConsole(COT_OTHER, "    %u because the term has an algebraic factor", algebraicCount);
}

void HyperCullenWoodallApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t x, y;
   int32_t  sign;
   uint64_t minPrime;
   
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

      if (sscanf(buffer, "%u %u %d", &x, &y, &sign) != 3)
         FatalError("Line %s is malformed", buffer);

      if (sign != +1 && sign != -1)
         FatalError("Only +1 or -1 is supported");

      if (haveBitMap)
      {
         if (sign == +1)
            iv_PlusTerms[BIT(x, y)] = true;
         else
            iv_MinusTerms[BIT(x, y)] = true;
      }
      else
      {
         if (sign == +1)
            ib_IsPlus = true;
         
         if (sign == -1)
            ib_IsMinus = true;
         
         if (ii_MinX == 0 || x < ii_MinX)
            ii_MinX = x;

         if (x > ii_MaxX)
            ii_MaxX = x;

         if (ii_MinY == 0 || y < ii_MinY)
            ii_MinY = y;

         if (y > ii_MaxY)
            ii_MaxY = y;
      }
         
      il_TermCount++;
   }

   fclose(fPtr);
}

bool HyperCullenWoodallApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint32_t x1, x2, y1, y2;
   int32_t  sign;
   
   if (sscanf(term, "%u^%u*%u^%u%d", &x1, &y1, &y2, &x2, &sign) != 5)
      FatalError("Could not parse term %s\n", term);

   if (x1 != x2)
      FatalError("x values for term %s do not match (%u != %u)\n", term, x1, x2);
   
   if (y1 != y2)
      FatalError("y values for term %s do not match (%u != %u)\n", term, y1, y2);

   if (sign != +1 && sign != -1)
      FatalError("sign term %s must be +1 or -1\n", term);
   
   if (x1 < ii_MinX || x1 > ii_MaxX)
      return false;
   
   if (y1 < ii_MinY || y1 > ii_MaxY)
      return false;

   if (!ib_IsMinus && sign != +1)
      FatalError("+/- for term %s is not the expected value\n", term);
   
   if (!ib_IsPlus && sign != -1)
      FatalError("+/- for term %s is not the expected value\n", term);

   VerifyFactor(theFactor, x1, y1, sign);

   uint64_t bit = BIT(x1, y1);
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
   uint32_t x, y;
   uint64_t bit, termsCounted = 0;
   
   ip_FactorAppLock->Lock();

   fPtr = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!fPtr)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   fprintf(fPtr, "ABC $a^$b*$b^$a$c // Sieved to %" PRIu64"\n", largestPrime);

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit])
         {
            fprintf(fPtr, "%u %u +1\n", x, y);
            termsCounted++;
         }
         
         if (iv_MinusTerms[bit])
         {
            fprintf(fPtr, "%u %u -1\n", x, y);
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
      snprintf(extraText, maxTextLength, "%d <= x <= %d, %d <= y <= %d (+1 only)",ii_MinX, ii_MaxX, ii_MinY, ii_MaxY);
   else if (!ib_IsPlus)
      snprintf(extraText, maxTextLength, "%d <= x <= %d, %d <= y <= %d (-1 only)",ii_MinX, ii_MaxX, ii_MinY, ii_MaxY);
   else
      snprintf(extraText, maxTextLength, "%d <= x <= %d, %d <= y <= %d (both +1 and -1)",ii_MinX, ii_MaxX, ii_MinY, ii_MaxY);
}

bool HyperCullenWoodallApp::ReportFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign)
{
   uint64_t bit;
   bool     removedTerm = false;
   
   if (x < ii_MinX || x > ii_MaxX)
      return false;
   
   if (y < ii_MinY || y > ii_MaxY)
      return false;
   
   if (sign != -1 && sign != +1)
      FatalError("expecting sign of +1 or -1");
   
   if (sign == +1 && !ib_IsPlus)
      return false;
   
   if (sign == -1 && !ib_IsMinus)
      return false;
   
   VerifyFactor(theFactor, x, y, sign);
   
   ip_FactorAppLock->Lock();

   bit = BIT(x, y);
   
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
      LogFactor(theFactor, "%u^%u*%u^%u%+d", x, y, y, x, sign);
   }
   
   ip_FactorAppLock->Release();

   return removedTerm;
}

void  HyperCullenWoodallApp::VerifyFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign)
{
   bool     isValid = true;
   
   MpArith  mp(theFactor);
   MpRes    resXY = mp.pow(mp.nToRes(x), y);
   MpRes    resYX = mp.pow(mp.nToRes(y), x);
   
   resXY = mp.mul(resXY, resYX);

   uint64_t res = mp.resToN(resXY);
   
   if (sign > 0 && res != theFactor - 1)
      isValid = false;
   
   if (sign < 0 && res != 1)
      isValid = false;
   
   if (isValid)
      return;
   
   FatalError("%" PRIu64" does not divide %u^%u*%u^%u%+d", theFactor, x, y, y, x, sign);
}

void  HyperCullenWoodallApp::BuildTerms(void)
{
   uint32_t bit;

   // This could be done using a lot less memory, but the time and effort to do that
   // is significant.  That will be left as a future task for anyone intested in doing it.
   ip_FactorAppLock->Lock();

   if (ip_xTerms == NULL)
      ip_xTerms = (base_t *) xmalloc(ii_MaxX - ii_MinX + 2, sizeof(base_t) , "xTerms");

   base_t *xTerm = ip_xTerms;

   for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
   {      
      xTerm->base = x;
      xTerm->powerCount = 0;
      xTerm->indexedByPower = false;

      for (uint32_t y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            xTerm->powerCount++;
      }

      xTerm->powers.resize(xTerm->powerCount);
      xTerm->residues.resize(xTerm->powerCount);
      
      xTerm->powerCount = 0;

      for (uint32_t y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            xTerm->powers[xTerm->powerCount] = y;
            xTerm->powerCount++;
         }
      }

      xTerm++;
   }

   if (ip_yTerms == NULL)
      ip_yTerms = (base_t *) xmalloc(ii_MaxY - ii_MinY + 2, sizeof(base_t) , "yTerms");
   
   base_t *yTerm = ip_yTerms;
   
   for (uint32_t y=ii_MinY; y<=ii_MaxY; y++)
   {
      yTerm->base = y;

      for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      {
         bit = BIT(x, y);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            yTerm->powerCount++;
      }

      yTerm->powers.resize(yTerm->powerCount);
      yTerm->residues.resize(ii_MaxX - ii_MinX + 1);
      yTerm->indexedByPower = true;

      yTerm->powerCount = 0;

      for (uint32_t x=ii_MinX; x<=ii_MaxX; x++)
      {
         bit = BIT(x, y);

         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            yTerm->powers[yTerm->powerCount] = x;
            yTerm->powerCount++;
         }
      }
      
      yTerm++;
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
//   x1 x1y1 x1y2 ... x1yn 0
//   x2 x2y1 x2y2 ... x2yn 0
//    ... ... ... ... ...
//   xm xmy1 xmy2 ... xmyn 0
//   0
// Each group of terms for x will end with a 0 indicating that the next term
// after the 0 is for the next x.  Two 0s in a row indicates end of the group.
gpugroup_t  *HyperCullenWoodallApp::GetGpuGroupsOfTerms(void)
{
   uint32_t    bit, x, y;
   uint32_t    termsForX = 0;
   gpugroup_t *fgPtr = NULL;
   gpugroup_t *cgPtr = NULL;
   uint32_t    index;
   uint32_t    distinctXY = 0;
   uint32_t   *terms;
   int         groups = 1;
   
   ip_FactorAppLock->Lock();

   index = 0;

   x = ii_MinX;
   
   fgPtr = (gpugroup_t *) xmalloc(1, sizeof(gpugroup_t), "terms");
   fgPtr->terms = (uint32_t *) xmalloc(ii_MaxGpuSteps, sizeof(uint32_t), "steps");
   fgPtr->nextGroup = NULL;

   cgPtr = fgPtr;

   terms = fgPtr->terms;
   do
   {
      termsForX = 0;
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            termsForX++;
      }

      if (termsForX == 0)
         continue;

      // If not enough space for all y for this x, then put this x into the next group.
      if ((termsForX + index) >= (ii_MaxGpuSteps - 2))
      {
         if (index == 0)
            FatalError("-S must be at least %u", termsForX+2);

         // Signal that there are no more x for this group
         terms[index] = 0;
         
         cgPtr->nextGroup = (gpugroup_t *) xmalloc(1, sizeof(gpugroup_t), "terms");
         
         cgPtr = (gpugroup_t *) cgPtr->nextGroup;
         cgPtr->terms = (uint32_t *) xmalloc(ii_MaxGpuSteps, sizeof(uint32_t), "steps");
         index = 0;
         
         terms = cgPtr->terms;
         groups++;
      }
      
      terms[index] = x;
      index++;

      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            terms[index] = y;
            index++;
            distinctXY++;
         }
      }

      // Signal that there are no more y for this x
      terms[index] = 0;
      index++;
      
      x++;
   } while (x <= ii_MaxX);

   ip_FactorAppLock->Release();

   WriteToConsole(COT_OTHER, "After building GPU groups.  %u groups are needed to hold %u distinct x/y for %llu terms", groups, distinctXY, il_TermCount);
   return fgPtr;
}
#endif
