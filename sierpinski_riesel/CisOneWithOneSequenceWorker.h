/* CisOneWithOneSequenceWorker.h -- (C) Mark Rodenkirch, January 2021

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CisOneWithOneSequenceWorker_H
#define _CisOneWithOneSequenceWorker_H

#include "SierpinskiRieselApp.h"
#include "CisOneWithOneSequenceHelper.h"
#include "AbstractWorker.h"
#include "../core/HashTable.h"
#include "../core/MpArith.h"

using namespace std;

class CisOneWithOneSequenceWorker : public AbstractWorker
{
public:
   CisOneWithOneSequenceWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper);

   ~CisOneWithOneSequenceWorker(void) {};

   void              Prepare(uint64_t largestPrimeTested, uint32_t bestQ);
   
   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);
   
private:   
   sp_t              GetParity(uint64_t p);
   
   void              TestSinglePrime(uint64_t p, sp_t parity);
   
   uint32_t          SetupDiscreteLog(MpArith mp, MpRes resBase, MpRes resInvBase, MpRes resNegCK, sp_t parity);
   
   void              BuildLookupsAndClimbLadder(MpArith mp, MpRes resBase, MpRes resNegCK, uint16_t cqIdx, uint16_t ssCount, uint16_t *seqQs);
   
   uint32_t          BabySteps(MpArith mp, MpRes resBase, MpRes resInvBase, uint32_t babySteps);

   CisOneWithOneSequenceHelper *ip_CisOneHelper;

   // See SierpinskiRieselApp.h to see how these are defined.
   uint32_t          ii_BaseMultiple;
   uint32_t          ii_LimitBase;
   uint32_t          ii_PowerResidueLcm;
   
   uint32_t          ii_SieveLow;
   
   HashTable        *ip_HashTable;
   
   MpRes            *resBJ;         // there is one per Q
   MpRes            *resBD;         // there is one per Q
   MpRes            *resX;
   
   MpRes             resBexpQ;
   
   int16_t          *ip_DivisorShifts;
   uint16_t         *ip_PowerResidueIndices;
   
   uint32_t         *ip_CongruentQIndices;
   uint32_t         *ip_LadderIndices;
   
   uint16_t         *ip_AllQs;
   uint16_t         *ip_AllLadders;
      
   legendre_t       *ip_Legendre;
   uint8_t          *ip_LegendreTable;
   
   uint32_t          ii_Dim1;
   uint32_t          ii_Dim2;
   uint32_t          ii_Dim3;
   uint32_t          ii_Dim4;
};

#endif

