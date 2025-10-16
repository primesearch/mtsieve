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
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint32_t          GetN(void) { return ii_N; };
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
   
   bool              ReportFactor(uint64_t theFactor, uint64_t k, bool verifyFactor);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void);
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return !ib_TestTerms; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   void              OuptutAdditionalConsoleMessagesUponFinish(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   void              VerifyFactor(uint64_t theFactor, uint64_t k);
   void              TestRemainingTerms(void);
   
   std::vector<bool> iv_MMPTerms0;
   std::vector<bool> iv_MMPTerms1;
   
   std::string       is_InputFileName;
   std::string       is_OutputFileName;

   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_N;
   
   bool              ib_TestTerms;
   uint32_t          ii_MaxGpuFactors;
};

#endif

