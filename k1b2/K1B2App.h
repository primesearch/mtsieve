/* AFSieveApp.h -- (C) Mark Rodenkirch, November 2016

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _K1B2App_H
#define _K1B2App_H

#include "../core/FactorApp.h"

class K1B2App : public FactorApp
{
public:
   K1B2App(void);

   ~K1B2App(void) {};

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };
   
   int64_t           GetMinC(void) { return il_MinC; };
   int64_t           GetMaxC(void) { return il_MaxC; };

   bool              ReportFactor(uint64_t theFactor, uint32_t n, int64_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested) {};
   bool              IsWritingOutputTermsFile(void){ return true; };
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   void              WriteOutputTermsFile(uint64_t largestPrime);
   uint64_t          WriteABCDTermsFile(char *fileName, uint32_t minN, uint64_t maxPrime);

private:  
   std::vector<std::vector<bool>>  iv_Terms;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   int64_t           il_MinC;
   int64_t           il_MaxC;
};

#endif

