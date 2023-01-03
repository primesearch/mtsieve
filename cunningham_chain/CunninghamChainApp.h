/* CunninghamChainApp.h -- (C) Mark Rodenkirch, June 2022

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CunninghamChainApp_H
#define _CunninghamChainApp_H

#include "../core/FactorApp.h"

#define KMAX_MAX (UINT64_C(1)<<62)
#define NMAX_MAX (1 << 31)

typedef enum { FF_UNKNOWN = 1, FF_CC, FF_NEWPGEN } format_t;
typedef enum { CCT_UNKNOWN = 0, CCT_FIRSTKIND, CCT_SECONDKIND } chaintype_t;
typedef enum { TT_UNKNOWN = 0, TT_BN, TT_PRIMORIAL, TT_FACTORIAL } termtype_t;

class CunninghamChainApp : public FactorApp
{
public:
   CunninghamChainApp(void);

   ~CunninghamChainApp(void);

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   chaintype_t       GetChainKind(void) { return it_ChainKind; };
   termtype_t        GetTermType(void) { return it_TermType; };
   uint32_t          GetChainLength(void) { return ii_ChainLength; };
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetN(void) { return ii_N; };
   uint64_t         *GetTerms(void) { return il_Terms; };
   
   void              ReportFactor(uint64_t theFactor, uint64_t k, uint32_t termInChain);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
 
   void              BuildPrimorialTerms(void);
   void              BuildFactorialTerms(void);
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   uint64_t          WriteCCTermsFile(uint64_t largestPrime, FILE *termsFile);
   uint64_t          WriteNewPGenTermsFile(uint64_t largestPrime, FILE *termsFile);
   
   void              VerifyFactor(uint64_t theFactor, uint64_t k, uint32_t termInChain);
   
   std::vector<bool> iv_Terms;
   
   format_t          it_Format;
   
   std::string       is_InputFileName;
   std::string       is_OutputFileName;

   chaintype_t       it_ChainKind;
   termtype_t        it_TermType;
   uint32_t          ii_ChainLength;
   
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
   
   uint32_t         *ii_Primes;
   uint64_t         *il_Terms;
};

#endif

