/* SierpinskiRieselApp.h -- (C) Mark Rodenkirch, October 2018

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SierpinskiRieselApp_H
#define _SierpinskiRieselApp_H

#include "../core/FactorApp.h"
#include "AbstractSubsequenceHelper.h"

#define NMAX_MAX (1 << 31)

typedef enum { FF_UNKNOWN = 1, FF_ABCD, FF_ABC, FF_BOINC, FF_NUMBER_PRIMES } format_t;

class SierpinskiRieselApp : public FactorApp
{
public:
   SierpinskiRieselApp(void);

   ~SierpinskiRieselApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t thePrime, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);

   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };
   uint64_t          GetSmallSievePrimeLimit(void) { return il_SmallPrimeSieveLimit; };
   bool              UseAvxIfAvailable(void) { return ib_UseAvx; };
      
   AbstractSubsequenceHelper   *GetAppHelper(void) { return ip_AppHelper; };
   
   uint32_t          GetSequenceCount(void) { return ii_SequenceCount; };
   
   void              AddSequence(uint64_t k, int64_t c, uint32_t d);
   
   seq_t            *GetSequence(uint64_t k, int64_t c, uint32_t d);

   void              ReportFactor(uint64_t thePrime, uint32_t seqIdx, uint32_t n);

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   seq_t            *ip_Sequences;
   AbstractSubsequenceHelper   *ip_AppHelper; 
   
   string            is_LegendreFileName;
   
   bool              LoadSequencesFromFile(char *fileName);
   void              ValidateAndAddNewSequence(char *arg);

   void              MakeSubsequences(bool newSieve);
   void              ComputeSmallSieveLimit(void);
      
   uint32_t          WriteABCDTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteABCTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteBoincTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteABCNumberPrimesTermsFile(seq_t *seq, uint64_t maxPrime, FILE *termsFile, bool allSequencesHaveDEqual1);
   
   bool              IsPrime(uint64_t p, uint64_t k, uint32_t n, int64_t c);
   bool              VerifyFactor(uint64_t thePrime, uint32_t seqIdx, uint32_t n);
   
   uint64_t          il_SmallPrimeSieveLimit;

   bool              ib_UseAvx;
   bool              ib_UseLengendreTables;
   bool              ib_HaveNewSequences;
   format_t          it_Format;
   
   uint32_t          ii_Base;
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;

   uint32_t          ii_SequenceCount;
   uint32_t          ii_SequenceCapacity;
};

#endif

