/* AFSieveApp.h -- (C) Mark Rodenkirch, September 2018

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _DMDivisorApp_H
#define _DMDivisorApp_H

#include <gmp.h>
#include "../core/FactorApp.h"

#define KMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

class DMDivisorApp : public FactorApp
{
public:
   DMDivisorApp(void);

   ~DMDivisorApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint32_t          GetN(void) { return ii_N; };
   
   bool              ReportFactor(uint64_t p, uint64_t k, bool verifyFactor);

protected:
   void              PreSieveHook(void);
   bool              PostSieveHook(void);
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   void              WriteOutputTermsFile(uint64_t largestPrime);
   bool              IsWritingOutputTermsFile(void){ return !ib_TestTerms; };
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   void              TestRemainingTerms(void);
   bool              IsDoubleMersenneDivisor(uint64_t k);
   void              CheckRedc(mp_limb_t *xp, uint32_t xn, uint32_t b, uint64_t k);
   bool              VerifyFactor(bool badFactorIsFatal, uint64_t thePrime, uint64_t k);
   
   vector<bool>      iv_MMPTerms;
   
   string            is_InputFileName;
   string            is_OutputFileName;

   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_N;
   
   bool              ib_TestTerms;
   uint64_t          il_KPerChunk;
   
   uint64_t          il_MinKOriginal;
   uint64_t          il_MaxKOriginal;
   uint64_t          il_MinKInChunk;
      
   uint64_t          il_StartSievingUS;
   uint64_t          il_TotalTerms;
   uint64_t          il_TotalTermsEvaluated;
   uint64_t          il_TotalTermsInChunk;
};

#endif

