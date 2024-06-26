/* SierpinskiRieselApp.h -- (C) Mark Rodenkirch, October 2018

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _SierpinskiRieselApp_H
#define _SierpinskiRieselApp_H

#include <set>

#include "../core/FactorApp.h"
#include "AbstractSequenceHelper.h"
#include "../core/BigHashTable.h"

#define NMAX_MAX (1 << 31)

typedef enum { FF_UNKNOWN = 1, FF_ABCD, FF_ABC, FF_BOINC, FF_NUMBER_PRIMES, FF_MFAKT } format_t;

class SierpinskiRieselApp : public FactorApp
{
public:
   SierpinskiRieselApp(void);

   ~SierpinskiRieselApp(void);

   void              Help(void);
   void              AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(uint64_t theFactor, const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText, uint32_t maxTextLength);

   uint32_t          GetBase(void) { return ii_Base; };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };
   uint64_t          GetSmallSievePrimeLimit(void) { return il_SmallPrimeSieveLimit; };
   uint32_t          GetSquareFreeBase(void) { return ii_SquareFreeB; };
   
   AbstractSequenceHelper   *GetAppHelper(void) { return ip_AppHelper; };
   
   void              AddSequence(uint64_t k, int64_t c, uint32_t d);
   
   seq_t            *GetSequence(uint64_t k, int64_t c, uint32_t d, seq_t *startSeqPtr);

   void              ReportFactor(uint64_t theFactor, seq_t *seqPtr, uint32_t n, bool verifyFactor);

   bool              CanUseCIsOneLogic(void) { return !ib_UseGenericLogic; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint64_t          HasSingleC(void) { return ib_HaveSingleC; };
   uint64_t          GetLegendreTableBytes(void) { return il_LegendreTableBytes; };
   std::string       GetLegendreDirectoryName(void) { return is_LegendreDirectoryName; };
   
   seq_t            *GetFirstSequenceAndSequenceCount(uint32_t &count) { count = ii_SequenceCount; return ip_FirstSequence; };
   uint32_t          GetSequenceCount(void) { return ii_SequenceCount; };
   
   double            GetGiantStepFactor(void) { return id_GiantStepFactor; };
   uint32_t          GetBaseMultipleMulitplier(void) { return ii_BaseMultipleMultiplier; };
   uint32_t          GetPowerResidueLcmMultiplier(void) { return ii_PowerResidueLcmMulitplier; };
   uint32_t          GetLimitBaseMultiplier(void) { return ii_LimitBaseMultiplier; };
   bool              ShowQEffort(void) { return ib_ShowQEffort; };
   uint32_t          GetUserBestQ(void) { return ii_UserBestQ; };
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          GetMaxGpuFactors(void) { return ii_MaxGpuFactors; };
   uint32_t          GetKernelCount(void) { return ii_KernelCount; };
   void              UseGpuWorkersUponRebuild(void) { ib_UseGPUWorkersUponRebuild = true; };
#endif

protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);
   void              OuptutAdditionalConsoleMessagesUponFinish(void);
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);

private:
   double            id_GiantStepFactor;
   uint32_t          ii_BaseMultipleMultiplier;
   uint32_t          ii_PowerResidueLcmMulitplier;
   uint32_t          ii_LimitBaseMultiplier;

   seq_t            *ip_FirstSequence;
   AbstractSequenceHelper   *ip_AppHelper; 
   
   uint64_t          il_LegendreTableBytes;
   std::string       is_LegendreDirectoryName;
   std::string       is_SequencesToRemove;
   
   std::set<uint64_t>     is_SequenceKs;
   std::set<std::string>  is_SequenceSet;
   seq_t            *ip_LastSequence;
   
   bool              LoadSequencesFromFile(char *fileName);
   void              ValidateAndAddNewSequence(char *arg);

   void              MakeSubsequences(bool newSieve, uint64_t largestPrimeTested);
   
   void              RemoveSequences(void);
   void              RemoveSequence(const char *sequence);
   void              RemoveSequence(uint64_t k, uint32_t b, int64_t c, uint32_t d);

   void              RemoveSequencesWithNoTerms(void);
   void              UpdateSequenceIndexes(void);
   
   void              CheckForLegendreSupport(void);
   void              RemoveN(void);
   void              CheckForUnsievableSequences(void);
      
   void              WriteOutputTermsFilesByQ(void);
   void              WriteSequenceFilesByQ(uint32_t q);
   uint32_t          WriteOutputTermsFile(uint64_t largestPrime, uint32_t q);
   uint32_t          WriteABCDTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteABCTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteBoincTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   uint32_t          WriteABCNumberPrimesTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile, bool allSequencesHaveDEqual1);
   uint32_t          WriteMfaktTermsFile(seq_t *seqPtr, uint64_t maxPrime, FILE *termsFile);
   
   bool              IsPrime(uint64_t p, seq_t *seqPtr, uint32_t n);
   void              VerifyFactor(uint64_t theFactor, seq_t *seqPtr, uint32_t n);
   
   uint64_t          GetSquareFreeFactor(uint64_t n, std::vector<uint64_t> primes);
      
#if defined(USE_OPENCL) || defined(USE_METAL)
   bool              ib_UseGPUWorkersUponRebuild;
#endif
   
   uint64_t          il_SmallPrimeSieveLimit;

   bool              ib_HaveNewSequences;
   bool              ib_HaveSingleC;
   bool              ib_HaveGenericWorkers;
   bool              ib_RemoveN;
   bool              ib_OnlyPrimeNs;
   bool              ib_UseGenericLogic;
   format_t          it_Format;
   double            id_EstimatedPrimes;
   
   uint32_t          ii_Base;
   uint32_t          ii_SquareFreeB;    // product of squery free factors of the base
   
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   uint64_t          il_MaxK;
   uint64_t          il_MaxAbsC;
   uint32_t          ii_MaxD;

   bool              ib_Algebraic;
   bool              ib_SplitByBestQ;
   bool              ib_ShowQEffort;
   uint32_t          ii_UserBestQ;
   uint32_t          ii_SequenceCount;
   
#if defined(USE_OPENCL) || defined(USE_METAL)
   uint32_t          ii_GpuFactorDensity;
   uint32_t          ii_MaxGpuFactors;
   uint32_t          ii_KernelCount;
#endif
};

#endif

