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
typedef enum { TT_UNKNOWN = 0, TT_BN, TT_PRIMORIAL, TT_FACTORIAL } termtype_t;

class TwinApp : public FactorApp
{
public:
   TwinApp(void);

   ~TwinApp(void);

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetN(void) { return ii_N; };
   termtype_t        GetTermType(void) { return it_TermType; };
   uint64_t         *GetTerms(void) { return il_Terms; };
   
   void              ReportFactor(uint64_t theFactor, uint64_t k, int32_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   
   void              BuildPrimorialTerms(void);
   void              BuildFactorialTerms(void);
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   void              OuptutAdditionalConsoleMessagesUponFinish(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   uint64_t          WriteABCDTermsFile(uint64_t maxPrime, FILE *termsFile);
   uint64_t          WriteABCTermsFile(uint64_t maxPrime, FILE *termsFile);
   uint64_t          WriteNewPGenTermsFile(uint64_t maxPrime, FILE *termsFile);
   void              AdjustMaxPrime(void);
   void              VerifyFactor(uint64_t theFactor, uint64_t k, int32_t c);

   std::vector<bool> iv_TwinTerms;
   std::vector<bool> iv_MinusTerms;
   std::vector<bool> iv_PlusTerms;
   termtype_t        it_TermType;
   
   std::string       is_InputFileName;
   std::string       is_OutputFileName;

   format_t          it_Format;
   bool              ib_HalfK;

   uint64_t          il_MaxPrimeForValidFactor;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
   bool              ib_OnlyTwins;
   bool              ib_Remove;
   
   uint32_t         *ii_Primes;
   uint64_t         *il_Terms;
};

#endif

