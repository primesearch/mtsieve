/* AFSieveApp.h -- (C) Mark Rodenkirch, September 2018

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _TwinApp_H
#define _TwinApp_H

#include "../core/FactorApp.h"

#define KMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

typedef enum { FF_UNKNOWN = 1, FF_ABCD, FF_ABC, FF_NEWPGEN } format_t;

class TwinApp : public FactorApp
{
public:
   TwinApp(void);

   ~TwinApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetN(void) { return ii_N; };
   
   bool              ReportFactor(uint64_t p, uint64_t k, int32_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   uint64_t          WriteABCDTermsFile(uint64_t maxPrime, FILE *termsFile);
   uint64_t          WriteABCTermsFile(uint64_t maxPrime, FILE *termsFile);
   uint64_t          WriteNewPGenTermsFile(uint64_t maxPrime, FILE *termsFile);
   void              AdjustMaxPrime(void);

   vector<bool>      iv_TwinTerms;
   vector<bool>      iv_MinusTerms;
   vector<bool>      iv_PlusTerms;
   
   string            is_InputFileName;
   string            is_OutputFileName;

   format_t          it_Format;
   uint64_t          il_MaxPrimeForValidFactor;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
   bool              ib_OnlyTwins;
   bool              ib_Remove;
};

#endif

