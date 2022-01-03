/* CisOneWithMultipleSequencesWorker.h -- (C) Mark Rodenkirch, January 2021

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CisOneWithMultipleSequencesWorker_H
#define _CisOneWithMultipleSequencesWorker_H

#include "SierpinskiRieselApp.h"
#include "CisOneWithMultipleSequencesHelper.h"
#include "AbstractWorker.h"
#include "../core/HashTable.h"
#include "../core/MpArith.h"

using namespace std;

typedef struct {
   seq_t            *seqPtr;
   subseq_t         *ssPtr;
   uint32_t          q;
   MpRes             resBDCK;
} useable_subseq_t;

class CisOneWithMultipleSequencesWorker : public AbstractWorker
{
public:
   CisOneWithMultipleSequencesWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper);

   ~CisOneWithMultipleSequencesWorker(void) {};

   void              Prepare(uint64_t largestPrimeTested, uint32_t bestQ);
   
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);
   
private:   
   sp_t              GetParity(uint64_t p);
   
   void              TestSinglePrime(uint64_t p);
   
   uint32_t          SetupDiscreteLog(MpArith mp, MpRes resBase, MpRes resInvBase, MpRes resNegCK);
   
   uint32_t          ClimbLadder(MpArith mp, MpRes resBase);
   uint32_t          SetupDiscreteLog(MpArith mp, MpRes resBase, MpRes resInvBase);
   
   uint32_t          GetShift0Subsequences(MpArith mp, uint64_t bm);
   uint32_t          GetShiftXSubsequences(MpArith mp, uint64_t pShift, uint64_t bm, uint32_t r);
   
   uint32_t          BabySteps(MpArith mp, MpRes resBase, MpRes resInvBase, uint32_t babySteps);

   CisOneWithMultipleSequencesHelper *ip_CisOneHelper;
 
   uint32_t          ii_SieveLow;
   
   HashTable        *ip_HashTable;
   
   MpRes            *resBJ;         // there is one per Q
   MpRes            *resBD;         // there is one per Q
   MpRes            *resX;
   
   MpRes             resBexpQ;

   // See SierpinskiRieselApp.h to see how these are defined.
   uint32_t          ii_BaseMultiple;
   uint32_t          ii_LimitBase;
   uint32_t          ii_PowerResidueLcm;   
   
   uint32_t          ii_MaxBabySteps;
   
   int16_t          *ip_DivisorShifts;
   uint16_t         *ip_PowerResidueIndices;
   uint16_t          ii_PowerResidueIndexCount;

   legendre_t       *ip_Legendre;
   uint8_t          *ip_LegendreTable;
   
   bool              ib_AllSequencesHaveLegendreTables;
   
   useable_subseq_t *ip_UsableSubsequences;
   
   uint32_t         *ip_CongruentSubseqIndices;
   
   uint32_t         *ip_AllSubseqs;
   uint16_t         *ip_AllLadders;

   uint32_t          ii_Dim1;
   uint32_t          ii_Dim2;
   uint32_t          ii_Dim3;
};

#endif

