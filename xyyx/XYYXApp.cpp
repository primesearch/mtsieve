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
#define APP_VERSION     "1.8"

#define BIT(x, y)       ((((x) - ii_MinX) * GetYCount()) + ((y) - ii_MinY))

// This is declared in App.h, but implemented here.  This means that App.h
// can remain unchanged if using the GPUSieve framework for other applications.
App *get_app(void)
{
   return new XYYXApp();
}

XYYXApp::XYYXApp(void) : FactorApp()
{
   SetBanner(APP_NAME " v" APP_VERSION ", a program to find factors numbers of the form x^y+y^x or x^y-y^x");
   SetLogFileName("xyyxsieve.log");

   ii_MinX = 0;
   ii_MaxX = 0;
   ii_MinY = 0;
   ii_MaxY = 0;
   ii_SplitYCount = 0;
   ii_SplitYValue = 0;
   ii_CpuWorkSize = 10000;
   ii_GpuSteps = 5000;
   ib_IsPlus = false;
   ib_IsMinus = false;
   SetAppMinPrime(3);
   ib_UseAvx = true;
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

   shortOpts += "x:X:y:Y:s:S:z:Z:D";

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
		
      case 'z':
         status = Parser::Parse(arg, 1, 1000000000, ii_SplitYCount);
         break;
         
      case 'Z':
         status = Parser::Parse(arg, 1, 1000000000, ii_SplitYValue);
         break;
         
      case 'D':
         ib_UseAvx = false;
         status = P_SUCCESS;
         break;
         
      case 's':
         char value;
         
         status = Parser::Parse(arg, "+-", value);
         
         ib_IsPlus = (value == '+');
         ib_IsMinus = (value == '-');
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
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);
      
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

      // Force x > y due to how searches are typically executed
      if (ib_IsPlus && ii_MaxY > ii_MaxX)
         FatalError("for + max x must be greater than max y");
      
      // Force x < y because we want x^y > y^x
      if (ib_IsMinus && ii_MaxY < ii_MaxX)
         FatalError("for - max x must be less than max y");

      il_TermCount = GetXCount() * GetYCount();
      
      iv_Terms.resize(il_TermCount);
      std::fill(iv_Terms.begin(), iv_Terms.end(), false);
         
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
   char     buffer[1000];
   uint32_t x, y;
   uint8_t  sign;
   uint64_t minPrime;
   
   if (!fPtr)
      FatalError("Unable to open input file %s", is_InputTermsFileName.c_str());

   if (fgets(buffer, sizeof(buffer), fPtr) == NULL)
      FatalError("File %s is empty", is_InputTermsFileName.c_str());

   if (sscanf(buffer, "ABC $a^$b%c$b^$a // Sieved to %" SCNu64"", &sign, &minPrime) != 2)
      FatalError("First line of the input file is malformed");

   ib_IsPlus = (sign == '+');
   ib_IsMinus = (sign == '-');
   
   if (!ib_IsMinus && !ib_IsPlus)
      FatalError("The sign must be '+' or '-'");
      
   // Reset this
   il_TermCount = 0;

   SetMinPrime(minPrime);

   while (fgets(buffer, sizeof(buffer), fPtr) != NULL)
   {
      if (!StripCRLF(buffer))
         continue;

      if (sscanf(buffer, "%u %u", &x, &y) != 2)
         FatalError("Line %s is malformed", buffer);

      if (y >= x)
         continue;

      if (haveBitMap)
      {
         iv_Terms[BIT(x, y)] = true;
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
      }
   }

   fclose(fPtr);
}

bool XYYXApp::ApplyFactor(uint64_t thePrime, const char *term)
{
   uint32_t x1, x2, y1, y2;
   uint8_t  c;
   
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
   if (iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_TermCount--;
      return true;
   }

   return false;
}

void XYYXApp::WriteOutputTermsFile(uint64_t largestPrime)
{
   FILE    *fPtr, *sPtr = NULL;
   uint32_t x, y, bit, yCount;
   uint64_t termsCounted = 0;
   char     termCountStr[50];
   char     termsCountedStr[50];


   ip_FactorAppLock->Lock();

   fPtr = fopen(is_OutputTermsFileName.c_str(), "w");

   if (!fPtr)
      FatalError("Unable to open input file %s", is_OutputTermsFileName.c_str());
   
   fprintf(fPtr, "ABC $a^$b%c$b^$a // Sieved to %" PRIu64"\n", (ib_IsPlus ? '+' : '-'), largestPrime);

   if (ii_SplitYCount > 0 || (ii_SplitYValue >= ii_MinY && ii_SplitYValue < ii_MaxY))
   {
      sPtr = fopen("xyyx_split.pfgw", "w");
      fprintf(sPtr, "ABC $a^$b%c$b^$a // Sieved to %" PRIu64"\n", (ib_IsPlus ? '+' : '-'), largestPrime);
   }
   
   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      yCount = 0;
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_Terms[bit])
            yCount++;
      }
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_Terms[bit])
         {
            if (yCount <= ii_SplitYCount || y <= ii_SplitYValue)
               fprintf(sPtr, "%u %u\n", x, y);
            else
               fprintf(fPtr, "%u %u\n", x, y);
            
            termsCounted++;
         }
      }
   }

   if (sPtr != NULL)
      fclose(sPtr);
      
   fclose(fPtr);
   
   if (termsCounted != il_TermCount)
   {
      sprintf(termCountStr, "%" PRIu64"", il_TermCount);
      sprintf(termsCountedStr, "%" PRIu64"", termsCounted);
      
      FatalError("Something is wrong.  Counted terms (%s) != expected terms (%s)", termsCountedStr, termCountStr);
   }
   
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
   
   if (ib_IsPlus && c == +1 && iv_Terms[bit])
   {
      iv_Terms[bit] = false;
      il_TermCount--;
      il_FactorCount++;
      removedTerm = true;
      LogFactor(p, "%u^%u+%u^%u", x, y, y, x);
   }
   
   if (ib_IsMinus && c == -1 && iv_Terms[bit])
   {
      iv_Terms[bit] = false;
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
   uint32_t   xGTEy = 0, yGTEx = 0;
   
   // Reset this
   il_TermCount = 0;

   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bool stillPlus = ib_IsPlus;
         bool stillMinus = ib_IsMinus;

         if (ib_IsMinus && x >= y)
         {
            xGTEy++;
            continue;
         }
         
         if (ib_IsPlus && y >= x)
         {
            yGTEx++;
            continue;
         }

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
   
         if (stillPlus || stillMinus)
         {
            iv_Terms[bit] = true;
            il_TermCount++;
         }
      }
   }

   WriteToConsole(COT_OTHER, "Quick elimination of terms info (in order of check):");
   
   if (ib_IsMinus)
      WriteToConsole(COT_OTHER, "    %u because x >= y", xGTEy);
   if (ib_IsPlus)
      WriteToConsole(COT_OTHER, "    %u because y >= x", yGTEx);
   
   WriteToConsole(COT_OTHER, "    %u because the term is even", evenCount);
   WriteToConsole(COT_OTHER, "    %u because x and y have a common divisor", commonDivisorCount);
}

// Build two tables for the worker.  One table will be used to hold x^y mod p.  One will be used to hold y^x mod p.
// This is a slightly more overhead to only having entries for remaining candidates, but it
// reduces memory usage by a factor of 10.
void   XYYXApp::GetTerms(uint32_t fpuRemaindersCount, uint32_t avxRemaindersCount, bases_t *bases)
{
   uint32_t    x;
   uint32_t    y;
   uint32_t    bit, idx, powerCount, powerIndex;
   uint64_t   *fpuRemaindersPtr = 0;
   double     *avxRemaindersPtr = 0;
   base_t     *xPowY, *yPowX;
   powerofx_t *powerOfX;
   powerofy_t *powerOfY;
   
   ip_FactorAppLock->Lock();

   bases->xPowY = xPowY = (base_t *) xmalloc((ii_MaxX - ii_MinX + 1) * sizeof(base_t));
   
   for (x=ii_MinX; x<=ii_MaxX; x++)
   {
      powerCount = 0;
      
      for (y=ii_MinY; y<=ii_MaxY; y++)
      {
         bit = BIT(x, y);
         
         if (iv_Terms[bit])
            powerCount++;
      }
      
      xPowY[x - ii_MinX].base = x;
      xPowY[x - ii_MinX].powerCount = powerCount;
         
      if (powerCount > 0)
      {
         xPowY[x - ii_MinX].powerIndices = (uint32_t *) xmalloc((ii_MaxY - ii_MinY + 1) * sizeof(uint32_t));
         xPowY[x - ii_MinX].powersOfX = (powerofx_t *) xmalloc(powerCount * sizeof(powerofx_t));
         
         // Allocate these only once for this x
         if (fpuRemaindersCount > 0)
            fpuRemaindersPtr = (uint64_t *) xmalloc(powerCount * fpuRemaindersCount * sizeof(uint64_t));
            
         if (avxRemaindersCount > 0)
            avxRemaindersPtr = (double *) xmalloc(powerCount * avxRemaindersCount * sizeof(double));

         xPowY[x - ii_MinX].fpuRemaindersPtr = fpuRemaindersPtr;
         xPowY[x - ii_MinX].avxRemaindersPtr = avxRemaindersPtr;
         
         powerIndex = 0;
         
         for (y=ii_MinY; y<=ii_MaxY; y++)
         {
            bit = BIT(x, y);
            
            if (iv_Terms[bit])
            {
               xPowY[x - ii_MinX].powerIndices[y - ii_MinY] = powerIndex;
               
               powerOfX = &xPowY[x - ii_MinX].powersOfX[powerIndex];
               
               powerOfX->y = y;
               
               if (fpuRemaindersCount)
               {
                  powerOfX->fpuRemainders = fpuRemaindersPtr;
                  fpuRemaindersPtr += fpuRemaindersCount;
               }
               
               if (avxRemaindersCount)
               {
                  powerOfX->avxRemainders = avxRemaindersPtr;
                  avxRemaindersPtr += avxRemaindersCount;
               }
               
               powerIndex++;
            }
         }
      }
   }  

   bases->yPowX = yPowX = (base_t *) xmalloc((ii_MaxY - ii_MinY + 1) * sizeof(base_t));
      
   for (y=ii_MinY; y<=ii_MaxY; y++)
   {
      powerCount = 0;
      
      for (x=ii_MinX; x<=ii_MaxX; x++)
      {
         bit = BIT(x, y);
         
         if (iv_Terms[bit])
            powerCount++;
      }
      
      yPowX[y - ii_MinY].base = y;
      yPowX[y - ii_MinY].powerCount = powerCount;
      
      if (powerCount > 0)
      {
         yPowX[y - ii_MinY].powersOfY = (powerofy_t *) xmalloc(powerCount * sizeof(powerofy_t));
         
         powerIndex = 0;
         
         for (x=ii_MinX; x<=ii_MaxX; x++)
         {
            bit = BIT(x, y);
            
            if (iv_Terms[bit])
            {
               powerOfY = (powerofy_t *) &yPowX[y - ii_MinY].powersOfY[powerIndex];
               
               powerOfY->x = x;
               
               // Link these together so that we can find x^y that corresponds to y^x
               idx = xPowY[x - ii_MinX].powerIndices[y - ii_MinY];
               
               powerOfY->powerOfX = &xPowY[x - ii_MinX].powersOfX[idx];

               powerIndex++;
            }            
         }
      }
   }
      
   for (x=ii_MinX; x<=ii_MaxX; x++)       
      if (xPowY[x - ii_MinX].powerCount > 0)
         xfree(xPowY[x - ii_MinX].powerIndices);
   
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
         
         if (iv_Terms[bit])
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
         
         if (iv_Terms[bit])
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
         
         if (iv_Terms[bit])
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
