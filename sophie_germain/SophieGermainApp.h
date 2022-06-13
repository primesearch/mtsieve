/* SophieGermainApp.h -- (C) Mark Rodenkirch, July 2020

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SophieGermainApp_H
#define _SophieGermainApp_H

#include "../core/FactorApp.h"

#define KMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

typedef enum { FF_UNKNOWN = 1, FF_ABCD, FF_NEWPGEN } format_t;

class SophieGermainApp : public FactorApp
{
public:
   SophieGermainApp(void);

   ~SophieGermainApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetN(void) { return ii_N; };
   bool              IsGeneralizedSearch(void) { return ib_GeneralizedSearch; };

   void              ReportFactor(uint64_t theFactor, uint64_t k, bool firstOfPair, bool verifyFactor);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   uint64_t          WriteABCDTermsFile(uint64_t maxPrime, FILE *termsFile);
   uint64_t          WriteNewPGenTermsFile(uint64_t maxPrime, FILE *termsFile);
   
   void              VerifyFactor(uint64_t theFactor, uint64_t k, bool firstOfPair);
   
   std::vector<bool> iv_Terms;
   
   bool              ib_GeneralizedSearch;
   format_t          it_Format;
   
   std::string       is_InputFileName;
   std::string       is_OutputFileName;

   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
};

#endif

