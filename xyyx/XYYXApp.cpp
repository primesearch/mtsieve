/* XYYXApp.cpp -- (C) Mark Rodenkirch, September 2012

   XYYXSieve/Wall-Sun-Sun Search OpenCL application

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <stdint.h>
#include <stdarg.h>
#include <time.h>
#include "../core/inline.h"
#include "XYYXApp.h"
#include "XYYXWorker.h"
#ifdef HAVE_GPU_WORKERS
#include "XYYXGpuWorker.h"
#endif

#define APP_NAME        "xyyxsieve"
#define APP_VERSION     "1.6"

#define BIT(x, y)       ((((x) - ii_MinX) * GetYCount()) + ((y) - ii_MinY))

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the GPUSieve framework for other applications.
App *get_app(void)
{
   return new XYYXApp();
}

XYYXApp::XYYXApp(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors numbers of the form x^y+y^x");
   SetLogFileName("xyyxsieve.log");

   ii_MinX = 0;
   ii_MaxX = 0;
   ii_MinY = 0;
   ii_MaxY = 0;
   ii_CpuWorkSize = 10000;
   ii_GpuSteps = 5000;
   ib_IsPlus = false;
   ib_IsMinus = false;
   SetAppMinPrime(3);
   ib_UseAvx = true;
   
  // SetBlockWhenProcessingFirstChunk(true);

#ifdef HAVE_GPU_WORKERS
   ib_SupportsGPU = true;
#endif
}

void XYYXApp::Help(void)
{
   FactorApp::ParentHelp();

   printf("-x --minx=x           minimum x to search\n");
   printf("-X --maxx=X           maximum x to search\n");
   printf("-y --miny=y           minimum y to search\n");
   printf("-Y --maxy=Y           maximum y to search\n");
   printf("-D --disableavx       disableavx\n");
   printf("-s --sign=+/-/b       sign to sieve for\n");
#ifdef HAVE_GPU_WORKERS
   printf("-S --step=S           steps iterated per call to GPU (default %d)\n", ii_GpuSteps);
#endif
}

void  XYYXApp::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
   FactorApp::ParentAddCommandLineOptions(shortOpts, longOpts);

   shortOpts += "x:X:y:Y:s:S:D";

   AppendLongOpt(longOpts, "minx",              required_argument, 0, 'x');
   AppendLongOpt(longOpts, "maxx",              required_argument, 0, 'X');
   AppendLongOpt(longOpts, "miny",              required_argument, 0, 'y');
   AppendLongOpt(longOpts, "sign",              required_argument, 0, 's');
   AppendLongOpt(longOpts, "disableavx",        no_argument, 0, 'D');
   AppendLongOpt(longOpts, "sign",              required_argument, 0, 's');
#ifdef HAVE_GPU_WORKERS
   AppendLongOpt(longOpts, "steps",             required_argument, 0, 'S');
#endif
}


parse_t XYYXApp::ParseOption(int opt, char *arg, const char *source)
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
		
      case 'D':
         ib_UseAvx = false;
         status = P_SUCCESS;
         break;
         
      case 's':
         char value;
         status = Parser::Parse(arg, "+-", value);
         if (value == '-')
            ib_IsMinus = true;
         if (value == '+')
            ib_IsPlus = true;
         break;
         
#ifdef HAVE_GPU_WORKERS
      case 'S':
         status = Parser::Parse(arg, 1, 1000000000, ii_GpuSteps);
         break;
#endif
   }

   return status;
}

void XYYXApp::ValidateOptions(void)
{
   if (is_OutputTermsFileName.length() == 0)
      is_OutputTermsFileName = "xyyx.pfgw";

   if (is_InputTermsFileName.length() > 0)
   {
      ProcessInputTermsFile(false);
      
      il_TermCount = GetXCount() * GetYCount();
      
      iv_PlusTerms.resize(il_TermCount);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
   
      iv_MinusTerms.resize(il_TermCount);
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

      if (ib_IsPlus && ii_MaxY > ii_MaxX)
         FatalError("for + max x must be greater than max y");
      
      if (ib_IsMinus && ii_MaxY < ii_MaxX)
         FatalError("for - max x must be less than max y");

      il_TermCount = GetXCount() * GetYCount();
      
      iv_PlusTerms.resize(il_TermCount);
      std::fill(iv_PlusTerms.begin(), iv_PlusTerms.end(), false);
   
      iv_MinusTerms.resize(il_TermCount);
      std::fill(iv_MinusTerms.begin(), iv_MinusTerms.end(), false);
      
      SetInitialTerms();
   }
   
   FactorApp::ParentValidateOptions();
}

Worker *XYYXApp::CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested)
{
   Worker *theWorker;

#ifdef HAVE_GPU_WORKERS  
   if (gpuWorker)
      theWorker = new XYYXGpuWorker(id, this);
   else
#endif
      theWorker = new XYYXWorker(id, this);
   
   return theWorker;
}

void XYYXApp::ProcessInputTermsFile(bool haveBitMap)
{
   FILE    *fPtr = fopen(is_InputTermsFileName.c_str(), "r");
   char     buffer[200];
   uint32_t x, y;
   int32_t  sign;
   uint64_t minPrime;
   
   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, 200,fPtr) == NULL)
      FatalError("File %s is empty", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC $a^$b$c*$b^$a // Sieved to %" SCNu64, &minPrime) != 1)
      FatalError("First line of the input file is malformed");

   // Reset this
   il_TermCount = 0;

   SetMinPrime(minPrime);

   while (fgets(buffer, 200, fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (sscanf(buffer, "%u %u %d", &x, &y, &sign) != 3)
         FatalError("Line %s is malformed", buffer);

      if (haveBitMap)
      {
         if (sign == +1)
            iv_PlusTerms[BIT(x, y)] = true;
         else
            iv_MinusTerms[BIT(x, y)] = true;
            
         il_TermCount++;
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

         if (sign == 1)
            ib_IsPlus = true;

         if (sign == -1)
            ib_IsMinus = true;
      }
   }

   fclose(fPtr);
}

bool XYYXApp::ApplyFactor(const char *term)
{
   uint32_t x1, x2, y1, y2;
   char     c;
   
   if (sscanf(term, "%u^%u%c%u^%u", &x1, &y1, &c, &y2, &x2) != 5)
      FatalError("Could not parse term %s\n", term);

   if (x1 != x2)
      FatalError("x values for term %s do not match (%u != %u)\n", term, x1, x2);
   
   if (y1 != y2)
      FatalError("y values for term %s do not match (%u != %u)\n", term, y1, y2);

   if (x1 < ii_MinX || x1 > ii_MaxX)
      return false;
   
   if (y1 < ii_MinY || y1 > ii_MaxY)
      return false;

   uint64_t bit = BIT(x1, y1);
   
   // No locking is needed because the Workers aren't running yet
   if (c == '+' && iv_PlusTerms[bit])
   {
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      return true;
   }
   
   if (c == '-' && iv_MinusTerms[bit])
   {
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      return true;
   }

   return false;
}

void XYYXApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *fPtr;
   uint32_t x, y, bit;
   uint64_t terms = 0;

   fPtr = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!fPtr)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   ip_FactorAppLock->Lock();

   fprintf(fPtr, "ABC $a^$b$c*$b^$a // Sieved to %" PRIu64"\n", largestPrime);

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (ib_IsPlus && iv_PlusTerms[bit])
         {
            fprintf(fPtr, "%u %u +1\n", x, y);
            terms++;
         }
         
         if (ib_IsMinus && iv_MinusTerms[bit])
         {
            fprintf(fPtr, "%u %u -1\n", x, y);
            terms++;
         }
      }
   }

   fclose(fPtr);
   
   if (terms != il_TermCount)
      FatalError("Something is wrong.  Counted terms (%u) != expected terms (%u)", terms, il_TermCount);
   
   ip_FactorAppLock->Release();
}

void  XYYXApp::GetExtraTextForSieveStartedMessage(char *extraText)
{   
   sprintf(extraText, "%d <= x <= %d, %d <= y <= %d",ii_MinX, ii_MaxX, ii_MinY, ii_MaxY);
}

bool XYYXApp::ReportFactor(uint64_t p, uint32_t x, uint32_t y, int32_t c)
{
   uint64_t bit;
   bool     removedTerm = false;
   
   if (x < ii_MinX || x > ii_MaxX)
      return false;
   
   if (y < ii_MinY || y > ii_MaxY)
      return false;
   
   ip_FactorAppLock->Lock();

   bit = BIT(x, y);
   
   if (ib_IsPlus && c == +1 && iv_PlusTerms[bit])
   {
      iv_PlusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      removedTerm = true;
      LogFactor(p, "%u^%u+%u^%u", x, y, y, x);
   }
   
   if (ib_IsMinus && c == -1 && iv_MinusTerms[bit])
   {
      iv_MinusTerms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      removedTerm = true;
      LogFactor(p, "%u^%u-%u^%u", x, y, y, x);
   }
   
   ip_FactorAppLock->Release();

   return removedTerm;
}


void  XYYXApp::SetInitialTerms(void)
{
   uint32_t   x, y, bit;
   uint32_t   evenCount = 0;
   uint32_t   commonDivisorCount = 0;
   
   // Reset this
   il_TermCount = 0;

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bool stillPlus = ib_IsPlus;
         bool stillMinus = ib_IsMinus;

         // If x and y are both odd, then both x^y+y^x and x^y-y^x are even
         // so we don't need to add them
         if (x & 1 && y & 1)
         {
            if (ib_IsPlus) evenCount++;
            if (ib_IsMinus) evenCount++;
            continue;
         }

         // If x and y are both even, then both x^y+y^x and x^y-y^x are even
         // so we don't need to add them
         if (!(x & 1) && !(y & 1))
         {
            if (ib_IsPlus) evenCount++;
            if (ib_IsMinus) evenCount++;
            continue;
         }
 
         // If they have a common divisor, then that common divisor will
         // divided x^y+y^x.
         if (gcd32(x, y) > 1)
         {
            if (ib_IsPlus) commonDivisorCount++;
            if (ib_IsMinus) commonDivisorCount++;
            continue;
         }

         bit = BIT(x, y);
   
         if (stillPlus)
         {
            iv_PlusTerms[bit] = true;
            il_TermCount++;
         }
         
         if (stillMinus)
         {
            iv_MinusTerms[bit] = true;
            il_TermCount++;
         }
      }
   }

   WriteToConsole(COT_OTHER, "Quick elimination of terms info (in order of check):");
   WriteToConsole(COT_OTHER, "    %u because the term is even", evenCount);
   WriteToConsole(COT_OTHER, "    %u because x and y have a common divisor", commonDivisorCount);
}

// Build two one dimensional arrays of terms.  The first array is will starts with the value for x
// and is followed by each y for that x that does not have a factor.  The list of y ends with a 0.
// Then the next x, and so on.  The list is done when it ends with two zeroes.  The second array
// is like the first, but it starts with the first y and has each x for that y.
//
// Imagine we have 7 x,y terms (100, 205), (100, 207), (100, 211), (100, 213), (101, 204), (101, 206), (102, 207)
// The xyTerms array would look like this: 100 205 207 211 213 0 101 204 206 0 102 207 0 0
// The yxTerms array would look like this: 204 101 0 206 101 0 207 100 102 0 211 100 0 213 100 0 0
void  XYYXApp::GetTerms(uint32_t *xyTerms, uint32_t *yxTerms)
{
   uint32_t x, y, bit;
   uint32_t *ptr;

   ip_FactorAppLock->Lock();
   
   ptr = xyTerms;
   
   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      *ptr = x;
      ptr++;
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit]) 
         {
            *ptr = y;
            ptr++;
         }
      }
      
      *ptr = 0;
      ptr++;
   }

   *ptr = 0;
   ptr++;
   
   ptr = yxTerms;
   
   for (y=ii_MinY; y<=ii_MaxY; y++)
   {
      *ptr = y;
      ptr++;
      
      for (x=ii_MinX; x<=ii_MaxX; x++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            *ptr = x;
            ptr++;
         }
      }
      
      *ptr = 0;
      ptr++;
   }

   *ptr = 0;
   ptr++;
   
   ip_FactorAppLock->Release();
   
}

#ifdef HAVE_GPU_WORKERS
uint32_t  XYYXApp::GetNumberOfGroups(void)
{
   uint32_t bit, x, y;
   uint32_t groupCount = 1;
   uint32_t termsInGroup = 0, termsForX = 0;
   bool     firstXInGroup;
   
   ip_FactorAppLock->Lock();

   firstXInGroup = true;
   termsInGroup = 0;

   x = ii_MinX;
   
   do
   {
      termsForX = 0;
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
         {
            termsInGroup++;
            termsForX++;
         }
      }
      
      // Account for x for the start of the group and 0 for the end of the group.
      termsInGroup += 2;
      
      // We need all y for the first x to fit into a group.
      if (firstXInGroup && termsForX >= (ii_GpuSteps - 2))
         FatalError("Too many terms for x = %u.  Increase setting for -S to %u and try again", x, termsForX+10);

      firstXInGroup = false;

      // If not enough space for all y for this x, then put this x into the next group.
      if ((termsForX + termsInGroup) >= (ii_GpuSteps - 2))
      {
         termsInGroup = 0;
         firstXInGroup = true;
         groupCount++;
         continue;
      }
      
      x++;
   } while (x <= ii_MaxX);

   ip_FactorAppLock->Release();
   
   return groupCount;
}

// Build multiple one dimensional arrays of terms in groups.
uint32_t   XYYXApp::GetGroupedTerms(uint32_t *terms)
{
   uint32_t  bit, x, y;
   uint32_t  groupIndex = 0;
   uint32_t  termsForX = 0;
   uint32_t  index, maxIndexForGroup;
   
   ip_FactorAppLock->Lock();

   index = 0;
   maxIndexForGroup = index + ii_GpuSteps;

   x = ii_MinX;
   
   do
   {
      termsForX = 0;
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_PlusTerms[bit] || iv_MinusTerms[bit])
            termsForX++;
      }

      // If not enough space for all y for this x, then put this x into the next group.
      if ((termsForX + index) >= (maxIndexForGroup - 2))
      {
         terms[index] = 0;
         groupIndex++;
         index = groupIndex * ii_GpuSteps;
         maxIndexForGroup = index + ii_GpuSteps;
         continue;
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
         }
      }

      terms[index] = 0;
      index++;
      
      x++;
   } while (x <= ii_MaxX);


   ip_FactorAppLock->Release();
   
   return groupIndex + 1;
}
#endif