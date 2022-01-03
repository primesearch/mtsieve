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
#include "AbstractSequenceHelper.h"

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
   uint32_t          GetSquareFreeBase(void) { return ii_SquareFreeB; };
   
   AbstractSequenceHelper   *GetAppHelper(void) { return ip_AppHelper; };
   
   void              AddSequence(uint64_t k, int64_t c, uint32_t d);
   
   seq_t            *GetSequence(uint64_t k, int64_t c, uint32_t d);

   void              ReportFactor(uint64_t thePrime, seq_t *seqPtr, uint32_t n, bool verifyFactor);

   bool              CanUseCIsOneLogic(void) { return ib_CanUseCIsOneLogic; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint64_t          GetLegendreTableBytes(void) { return il_LegendreTableBytes; };
   string            GetLegendreDirectoryName(void) { return is_LegendreDirectoryName; };
   
   seq_t            *GetFirstSequenceAndSequenceCount(uint32_t &count) { count = ii_SequenceCount; return ip_FirstSequence; };
   
   uint32_t          GetBaseMultipleMulitplier(void) { return ii_BaseMultipleMultiplier; };
   uint32_t          GetPowerResidueLcmMultiplier(void) { return ii_PowerResidueLcmMulitplier; };
   uint32_t          GetLimitBaseMultiplier(void) { return ii_LimitBaseMultiplier; };
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
   uint32_t          GetSequecesPerKernel(void) { return ii_SequencesPerKernel; };
   void              UseGpuWorkersUponRebuild(void) { ib_UseGPUWorkersUponRebuild = true; };
#endif

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   uint32_t          ii_BaseMultipleMultiplier;
   uint32_t          ii_PowerResidueLcmMulitplier;
   uint32_t          ii_LimitBaseMultiplier;

   seq_t            *ip_FirstSequence;
   AbstractSequenceHelper   *ip_AppHelper; 
   
   uint64_t          il_LegendreTableBytes;
   string            is_LegendreDirectoryName;
   string            is_SequencesToRemove;
   
   bool              LoadSequencesFromFile(char *fileName);
   void              ValidateAndAddNewSequence(char *arg);

   void              MakeSubsequences(bool newSieve, uint64_t largestPrimeTested);
   
   void              RemoveSequences(void);
   void              RemoveSequence(const char *sequence);
   void              RemoveSequence(uint64_t k, uint32_t b, int64_t c, uint32_t d);
   
   void              RemoveSequencesWithNoTerms(void);
   void              CheckForLegendreSupport(void);
      
   uint32_t          WriteABCDTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteABCTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteBoincTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteABCNumberPrimesTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile, bool allSequencesHaveDEqual1);
   
   bool              IsPrime(uint64_t p, seq_t *seqPtr, uint32_t n);
   bool              VerifyFactor(bool badFactorIsFatal, uint64_t thePrime, seq_t *seqPtr, uint32_t n);
   
   uint64_t          GetSquareFreeFactor(uint64_t n, vector<uint64_t> primes);
   
   bool              ib_CanUseCIsOneLogic;
   
#ifdef HAVE_GPU_WORKERS
   bool              ib_UseGPUWorkersUponRebuild;
#endif
   
   uint64_t          il_SmallPrimeSieveLimit;

   bool              ib_HaveNewSequences;
   format_t          it_Format;
   
   uint32_t          ii_Base;
   uint32_t          ii_SquareFreeB;    // product of squery free factors of the base
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   uint64_t          il_MaxK;
   uint64_t          il_MaxAbsC;

   uint32_t          ii_SequenceCount;
   
#ifdef HAVE_GPU_WORKERS
   uint32_t          ii_GpuFactorDensity;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_SequencesPerKernel;
#endif
};

#endif

