/* AFSieveApp.h -- (C) Mark Rodenkirch, November 2016

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _KBBApp_H
#define _KBBApp_H

#include "../core/AlgebraicFactorApp.h"

#define BMAX_MAX (UINT64_C(1)<<31)

#define BIT(b)          ((b) - ii_MinB)

class KBBApp : public AlgebraicFactorApp
{
public:
   KBBApp(void);

   ~KBBApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint64_t          GetK(void) { return il_K; };
   uint32_t          GetMinB(void) { return ii_MinB; };
   uint32_t          GetMaxB(void) { return ii_MaxB; };
   void              GetBases(uint32_t *bases);
   
   bool              ReportFactor(uint64_t prime, uint32_t b, int32_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   void              EliminateGfnAndMersenneTerms(void);
   void              EliminateTermsWithAlgebraicFactors(void);
   uint32_t          EliminateTermsWithSimpleRoots(void);

#ifdef __MINGW_PRINTF_FORMAT
   bool              CheckAlgebraicFactor(uint32_t b, int32_t c, const char *fmt, ...)  __attribute__ ((format (__MINGW_PRINTF_FORMAT, 4, 5)));
#else
   bool              CheckAlgebraicFactor(uint32_t b, int32_t c, const char *fmt, ...)  __attribute__ ((format (printf, 4, 5)));
#endif

   vector<bool>      iv_PlusTerms;
   vector<bool>      iv_MinusTerms;

   string            is_InputFileName;
   string            is_OutputFileName;

   uint64_t          il_K;
   uint32_t          ii_MinB;
   uint32_t          ii_MaxB;
};

#endif

