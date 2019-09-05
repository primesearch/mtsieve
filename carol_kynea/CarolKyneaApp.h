/* CarolKyneaApp.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CarolKyneaApp_H
#define _CarolKyneaApp_H

#include "../core/FactorApp.h"
#include "../core/SharedMemoryItem.h"

class CarolKyneaApp : public FactorApp
{
public:
   CarolKyneaApp(void);

   ~CarolKyneaApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);

   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };

   bool              ReportFactor(uint64_t p, uint32_t n, int32_t c);
   void              ReportPrime(uint64_t p, uint32_t n, int32_t c);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   
private:
   vector<bool>      iv_PlusTerms;
   vector<bool>      iv_MinusTerms;

   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
};

#endif

