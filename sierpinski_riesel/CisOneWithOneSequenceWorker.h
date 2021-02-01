/* CisOneWithOneSequenceWorker.h -- (C) Mark Rodenkirch, January 2021

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CisOneWithOneSequenceWorker_H
#define _CisOneWithOneSequenceWorker_H

#include "SierpinskiRieselApp.h"
#include "CisOneSequenceHelper.h"
#include "AbstractWorker.h"
#include "../core/HashTable.h"
#include "../core/MpArith.h"

#define L_BYTES(x) (((1+x)>>3)+1)
#define L_BYTE(x)  ((x)>>3)
#define L_BIT(x)   (1<<((x)&7))

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
   
   uint32_t          SetupDiscreteLog(MpArith mp, MpRes resBase, uint64_t resInvBase, sp_t parity);
   
   uint32_t          BuildHashTableAndClimbLadder(MpArith mp, MpRes resBase, MpRes resNegCK);
   
   uint32_t          BabySteps(MpArith mp, MpRes resBase, MpRes resInvBase, uint32_t babySteps);

   CisOneSequenceHelper *ip_CisOneHelper;
      
   uint32_t          ii_SieveLow;
   
   legendre_t       *ip_Legendre;
   HashTable        *ip_HashTable;
   
   MpRes            *resBJ;         // there is one per Q
   MpRes            *resBD;         // there is one per Q
   MpRes            *resX;
   
   MpRes             resBexpQ;
   
   bool              ib_UseLegendreTables;

   int16_t          *ip_DivisorShifts;
   uint16_t         *ip_PrlIndices;
   uint32_t          ii_PrlCount;
   
   uint16_t         *ip_Qs;
   uint16_t         *ip_Ladders;
};

#endif

