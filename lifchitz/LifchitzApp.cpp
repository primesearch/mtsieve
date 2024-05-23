/* LifchitzApp.cpp -- (C) Mark Rodenkirch, April 2024

   Sieve x^x + y^y or x^x-y^y

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

#include "LifchitzApp.h"
#include "LifchitzWorker.h"

#if defined(USE_OPENCL) || defined(USE_METAL)
#include "LifchitzGpuWorker.h"
#endif

#define APP_NAME        "lifsieve"
#define APP_VERSION     "1.7.1"

int sortByXY(const void *a, const void *b)
{
   term_t *aPtr = (term_t *) a;
   term_t *bPtr = (term_t *) b;

   if (aPtr->x == bPtr->x)
      return (aPtr->y - bPtr->y);
   
   return (aPtr->x - bPtr->x);
}

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the GPUSieve framework for other applications.
App *get_app(void)
{
   return new LifchitzApp();
}

LifchitzApp::LifchitzApp(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors numbers of the form x^x+y^y and x^x-y^y");
   SetLogFileName("lifsieve.log");

   ii_MinX = 0;
   ii_MaxX = 0;
   ii_MinY = 0;
   ii_MaxY = 0;
   ii_CpuWorkSize = 10000;
   ib_IsPlus = false;
   ib_IsMinus = false;
   SetAppMinPrime(3);
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   ii_XPerChunk = 1000000000;
   ii_MaxGpuFactors = GetGpuWorkGroups() * 10000;
#endif
}

void LifchitzApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-x --minx=x           minimum x to search\n");
   printf("-X --maxx=X           maximum x to search\n");
   printf("-y --miny=y           minimum y to search\n");
   printf("-Y --maxy=Y           maximum y to search\n");
   printf("-s --sign=+/-/b       sign to sieve for\n");
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   printf("-z --xperchunk=z      number of x per chunk per call to GPU (default %d)\n", ii_XPerChunk);

   printf("-M --maxfactors=M     max number of factors to support per GPU worker chunk (default %u)\n", ii_MaxGpuFactors);
#endif
}

void  LifchitzApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "x:X:y:Y:s:z:M:";

   AppendLongOpt(longOpts, "minx",              required_argument, 0, 'x');
   AppendLongOpt(longOpts, "maxx",              required_argument, 0, 'X');
   AppendLongOpt(longOpts, "miny",              required_argument, 0, 'y');
   AppendLongOpt(longOpts, "maxy",              required_argument, 0, 'Y');
   AppendLongOpt(longOpts, "sign",              required_argument, 0, 's');
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   AppendLongOpt(longOpts, "xchunks",           required_argument, 0, 'z');
   AppendLongOpt(longOpts, "maxfactors",        required_argument, 0, 'M');
#endif
}


parse_t LifchitzApp::ParseOption(int opt, char *arg, const char *source)
{
   parse_t status = P_UNSUPPORTED;

   status = FactorApp::ParentParseOption(opt, arg, source);
   if (status != P_UNSUPPORTED) return status;

   switch (opt)
   {
      case 'x':
         status = Parser::Parse(arg, 20, 1000000000, ii_MinX);
         break;

      case 'X':
         status = Parser::Parse(arg, 20, 1000000000, ii_MaxX);
         break;
		 
      case 'y':
         status = Parser::Parse(arg, 2, 1000000000, ii_MinY);
         break;

      case 'Y':
         status = Parser::Parse(arg, 2, 1000000000, ii_MaxY);
         break;
		
      case 's':
         char value;
         
         status = Parser::Parse(arg, "+-b", value);
         
         if (value == '+' || value == 'b')
            ib_IsPlus = true;
         
         if (value == '-' || value == 'b')
            ib_IsMinus = true;
         break;
         
#if defined(USE_OPENCL) || defined(USE_METAL)
      case 'z':
         status = Parser::Parse(arg, 1, 1000000000, ii_XPerChunk);
         break;

      case 'M':
         status = Parser::Parse(arg, 10, 1000000000, ii_MaxGpuFactors);
         break;
#endif
   }

   return status;
}

void LifchitzApp::ValidateOptions(void)
{
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "lifchitz.pfgw";

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      // This will be resized when after reading the file as multiple rows of
      // input can be collapsed into a single entry in the array.
      ip_Terms = (term_t *) xmalloc(il_TermCount, sizeof(term_t), "terms");
      
      il_TermCount = 0;
      
      ProcessInputTermsFile(true);
   }
   else
   {
      if (ii_MinX == 0)
         ii_MinX = 20;
	  
      if (ii_MaxX == 0)
         FatalError("max x has not been specified");

      if (ii_MinY == 0)
         ii_MinY = 2;

      if (ii_MaxY == 0)
         ii_MaxY = ii_MaxX - 1;

      if (ii_MinX > ii_MaxX)
         FatalError("min x must be less than or equal to max x");

      if (ii_MinY > ii_MaxY)
         FatalError("min y must be less than or equal to max y");

      // Force x > y due to how searches are typically executed
      if (ii_MaxY > ii_MaxX)
         FatalError("max x must be greater than max y");

      if (!ib_IsMinus && !ib_IsPlus)
         ib_IsMinus = ib_IsPlus = true;
      
      // Half of the terms are even, so we don't need to allocate memory for them.
      il_TermCount = 10 + (ii_MaxX - ii_MinX + 1) * (ii_MaxY - ii_MinY + 1);
      
      ip_Terms = (term_t *) xmalloc(il_TermCount, sizeof(term_t), "terms");
      
      SetInitialTerms();
   }
   
   CollapseTerms();

   FactorApp::ParentValidateOptions();

   // This will sieve beyond the limit, but we want to make sure that at least one prime
   // larger than this limit is passed to the worker even if the worker does not test it.
   SetMaxPrimeForSingleWorker(1000);
}

Worker *LifchitzApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#if defined(USE_OPENCL) || defined(USE_METAL)  
   if (gpuWorker)
      theWorker = new LifchitzGpuWorker(id, this);
   else
#endif
      theWorker = new LifchitzWorker(id, this);
   
   return theWorker;
}

void LifchitzApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[1000];
   uint32_t x, y;
   int32_t  sign;
   uint64_t minPrime, taIdx = 0;
   
   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("File %s is empty", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC $a^$a$b*$c^$c // Sieved to %" SCNu64"", &minPrime) != 1)
      FatalError("First line of the input file is malformed");
     
   // Reset this
   il_TermCount = 0;

   SetMinPrime(minPrime);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (!memcmp(buffer, "ABC $a^$a$b*$c^$c //", 20))
         continue;

      if (sscanf(buffer, "%u %d %u", &x, &sign, &y) != 3)
         FatalError("Line %s is malformed", buffer);

      if (sign != 1 && sign != -1)
         FatalError("Expected +1 or -1 on line %s", buffer);

      if (sign == +1) ib_IsPlus = true;
      if (sign == -1) ib_IsMinus = true;
         
      if (haveBitMap)
      {
         ip_Terms[taIdx].x = x;
         ip_Terms[taIdx].y = y;
         ip_Terms[taIdx].signs = 0;
         if (sign == -1) ip_Terms[taIdx].signs |= M_ONE;
         if (sign == +1) ip_Terms[taIdx].signs |= P_ONE;
         taIdx++;
      }
      else
      {
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

bool LifchitzApp::ApplyFactor(uint64_t theFactor, const char *term)
{
   uint32_t x1, x2, y1, y2;
   char     c;
   int32_t  sign;
   
   if (sscanf(term, "%u^%u%c%u^%u", &x1, &x2, &c, &y1, &y2) != 5)
      FatalError("Could not parse term %s\n", term);

   if (x1 != x2)
      FatalError("x values for term %s do not match (%u != %u)\n", term, x1, x2);
   
   if (y1 != y2)
      FatalError("y values for term %s do not match (%u != %u)\n", term, y1, y2);

   if (x1 < ii_MinX || x1 > ii_MaxX)
      return false;
   
   if (y1 < ii_MinY || y1 > ii_MaxY)
      return false;

   if (c != '+' && c != '-')
      FatalError("+/- for term %s is not the expected value\n", term);
   
   sign = (c == '+' ? +1 : -1);
   
   VerifyFactor(theFactor, x1, y1, sign);
   
   term_t key;
   key.x = x1;
   key.y = y1;

   term_t *entry = (term_t *) bsearch(&key, ip_Terms, il_TermArraySize, sizeof(term_t), sortByXY);
   
   return RemoveTerm(theFactor, entry, sign, false);
}

void LifchitzApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *fPtr;
   uint64_t termsCounted = 0;
   
   ip_FactorAppLock->Lock();

   fPtr = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!fPtr)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   fprintf(fPtr, "ABC $a^$a$b*$c^$c // Sieved to %" PRIu64"\n", largestPrime);

   for (uint32_t idx=0; idx<il_TermArraySize; idx++)
   {
      if (ip_Terms[idx].signs & M_ONE)
      {
         fprintf(fPtr, "%u -1 %u\n", ip_Terms[idx].x, ip_Terms[idx].y);
         termsCounted++;
      }

      if (ip_Terms[idx].signs & P_ONE)
      {
         fprintf(fPtr, "%u +1 %u\n", ip_Terms[idx].x, ip_Terms[idx].y);
         termsCounted++;
      }
   }

   fclose(fPtr);
   
   if (termsCounted != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", termsCounted, il_TermCount);
   
   ip_FactorAppLock->Release();
}

void  LifchitzApp::GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength)
{   
   snprintf(extraText, maxTextLength, "%u <= x <= %u, %u <= y <= %u, sign %s",ii_MinX, ii_MaxX, ii_MinY, ii_MaxY,
      ((ib_IsPlus & ib_IsMinus) ? "+/-" : (ib_IsPlus ? "+" : "-")));
}

bool LifchitzApp::ReportFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign)
{
   bool     removedTerm;
   bool     needToLock = (theFactor > GetMaxPrimeForSingleWorker());
   term_t  *entry;
   
   if (x < ii_MinX || x > ii_MaxX)
      return false;
   
   if (y < ii_MinY || y > ii_MaxY)
      return false;
      
   VerifyFactor(theFactor, x, y, sign);
   
   if (needToLock)
      ip_FactorAppLock->Lock();
   
   term_t key;
   key.x = x;
   key.y = y;

   entry = (term_t *) bsearch(&key, ip_Terms, il_TermArraySize, sizeof(term_t), sortByXY);

   removedTerm = RemoveTerm(theFactor, entry, sign, true);
      
   if (needToLock)
      ip_FactorAppLock->Release();

   return removedTerm;
}

bool LifchitzApp::RemoveTerm(uint64_t theFactor, term_t *term, int32_t sign, bool logFactor)
{
   bool needFactor = false;

   if (!term)
      return false;

   if (sign == -1 && term->signs & M_ONE) needFactor = true;
   if (sign == +1 && term->signs & P_ONE) needFactor = true;

   if (needFactor)
   {
      if (sign == -1) term->signs &= ~M_ONE;
      if (sign == +1) term->signs &= ~P_ONE;

      il_TermCount--;
      il_FactorCount++;
      
      if (logFactor)
         LogFactor(theFactor, "%u^%u%c%u^%u", term->x, term->x, ((sign == +1) ? '+' : '-'), term->y, term->y);
   }
   
   return needFactor;
}

void  LifchitzApp::VerifyFactor(uint64_t theFactor, uint32_t x, uint32_t y, int32_t sign)
{
   uint64_t xPowX, yPowY;
   bool     isValid = true;
   
   MpArith  mp(theFactor);
   MpRes    resXX = mp.pow(mp.nToRes(x), x);
   MpRes    resYY = mp.pow(mp.nToRes(y), y);
   
   xPowX = mp.resToN(resXX);
   yPowY = mp.resToN(resYY);
   
   if (sign == -1 && xPowX != yPowY)
      isValid = false;
   
   if (sign == +1 && xPowX + yPowY != theFactor)
      isValid = false;
   
   if (isValid)
      return;
   
   FatalError("%" PRIu64" does not divide %u^%u%c%u^%u", theFactor, x, x, ((sign == 1) ? '+' : '-'), y, y);
}

void  LifchitzApp::SetInitialTerms(void)
{
   uint32_t   x, y;
   uint32_t   evenCount = 0;
   uint32_t   commonDivisorCount = 0;
   uint32_t   yGTEx = 0;
   
   // Reset this
   il_TermArraySize = 0;
   il_TermCount = 0;

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         if (y >= x)
         {
            if (ib_IsPlus) yGTEx++;
            if (ib_IsMinus) yGTEx++;
            continue;
         }

         // If x and y are both odd, then both x^x+y^y and x^x-y^y are even
         // so we don't need to add them
         if (x & 1 && y & 1)
         {
            if (ib_IsPlus) evenCount++;
            if (ib_IsMinus) evenCount++;
            continue;
         }

         // If x and y are both even, then both x^x+y^y and x^x-y^y are even
         // so we don't need to add them
         if (!(x & 1) && !(y & 1))
         {
            if (ib_IsPlus) evenCount++;
            if (ib_IsMinus) evenCount++;
            continue;
         }
 
         // If they have a common divisor, then that common divisor will divide the term
         if (gcd32(x, y) > 1)
         {
            if (ib_IsPlus) commonDivisorCount++;
            if (ib_IsMinus) commonDivisorCount++;
            continue;
         }

         if (ib_IsMinus || ib_IsPlus)
         {
            ip_Terms[il_TermArraySize].x = x;
            ip_Terms[il_TermArraySize].y = y;
            ip_Terms[il_TermArraySize].signs = 0;

            if (ib_IsMinus)
            {
               ip_Terms[il_TermArraySize].signs |= M_ONE;
               il_TermCount++;
            }
            
            if (ib_IsPlus)
            {
               ip_Terms[il_TermArraySize].signs |= P_ONE;
               il_TermCount++;
            }
               
            il_TermArraySize++;
         }
      }
   }

   WriteToConsole(COT_OTHER, "Quick elimination of terms info (in order of check):");
   
   WriteToConsole(COT_OTHER, "    %u because y >= x", yGTEx);
   WriteToConsole(COT_OTHER, "    %u because the term is even", evenCount);
   WriteToConsole(COT_OTHER, "    %u because x and y have a common divisor", commonDivisorCount);
}

void      LifchitzApp::CollapseTerms(void)
{
   uint64_t  newTermArraySize = 0, taIdx, newTermCount;
   uint32_t  prevX = 0, prevY = 0;
   term_t   *newTerms;

   // sort by ascending x then ascending x then ascending y
   qsort(ip_Terms, il_TermCount, sizeof(term_t), sortByXY);

   for (uint64_t idx=0; idx<il_TermCount; idx++)
   {
      if (prevX == ip_Terms[idx].x && prevY == ip_Terms[idx].y)
         continue;

      newTermArraySize++;
      prevX = ip_Terms[idx].x;
      prevY = ip_Terms[idx].y;
   }
   
   newTerms = (term_t *) xmalloc(1 + newTermArraySize, sizeof(term_t), "terms");

   prevX = prevY = 0;
   newTermCount = taIdx = 0;
   
   for (uint64_t idx=0; idx<il_TermCount; idx++)
   {
      if (prevX == ip_Terms[idx].x && prevY == ip_Terms[idx].y)
      {
         if (ip_Terms[idx].signs & M_ONE)
         {
            newTermCount++;
            newTerms[taIdx-1].signs |= M_ONE;
         }

         if (ip_Terms[idx].signs & P_ONE)
         {
            newTermCount++;
            newTerms[taIdx-1].signs |= P_ONE;
         }
         
         continue;
      }

      newTermCount++;
      newTerms[taIdx].x = ip_Terms[idx].x;
      newTerms[taIdx].y = ip_Terms[idx].y;
      newTerms[taIdx].signs = ip_Terms[idx].signs;
      
      taIdx++;
      
      prevX = ip_Terms[idx].x;
      prevY = ip_Terms[idx].y;
   }
   
   if (il_TermCount != newTermCount)
      FatalError("Error encountered collapsing input terms.  Counted terms (%" PRIu64") != expected terms (%" PRIu64")", newTermCount, il_TermCount);

   xfree(ip_Terms);
   ip_Terms = newTerms;
   il_TermArraySize = newTermArraySize;
}
