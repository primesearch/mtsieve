/* AFSieveApp.h -- (C) Mark Rodenkirch, November 2016

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _FixedKBNApp_H
#define _FixedKBNApp_H

#include "../core/FactorApp.h"

#define CMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

#define BIT(c)          ((c) - il_MinC)

class FixedKBNApp : public FactorApp
{
public:
   FixedKBNApp(void);

   ~FixedKBNApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint64_t          GetK(void) { return il_K; };
   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetN(void) { return ii_N; };
   int64_t           GetMinC(void) { return il_MinC; };
   int64_t           GetMaxC(void) { return il_MaxC; };
   
   bool              ReportFactor(uint64_t p, int64_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void) {};
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:   
   vector<bool>      iv_Terms;

   string            is_Sequence;
   string            is_InputFileName;
   string            is_OutputFileName;

   uint64_t          il_K;
   int64_t           il_MinC;
   int64_t           il_MaxC;
   uint32_t          ii_Base;
   uint32_t          ii_N;
};

#endif

