/* CunninghamChainWorker.h -- (C) Mark Rodenkirch, June 2022

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CunninghamChainWorker_H
#define _CunninghamChainWorker_H

#include "CunninghamChainApp.h"
#include "../core/Worker.h"

using namespace std;

class CunninghamChainWorker : public Worker
{
public:
   CunninghamChainWorker(uint32_t myId, App *theApp);

   ~CunninghamChainWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   void              NotifyPrimeListAllocated(uint32_t primesInList);

private:
   void              TestSmallB(void);
   
   void              RemoveTerms(uint64_t thePrime, uint64_t k, uint32_t termInChain);

   void              BuildBaseInverses(void);
   uint32_t          EuclidExtendedGCD(uint32_t a, uint32_t base);
   
   CunninghamChainApp *ip_CunninghamChainApp;

   uint32_t         *ii_BaseInverses;
   uint64_t         *il_MyPrimeList;
   uint32_t         *ii_InverseList;
   
   chainkind_t       it_ChainKind;
   termtype_t        it_TermType;
   uint32_t          ii_ChainLength;
   
   bool              ib_HalfK;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_Base;
   uint32_t          ii_N;
   
   uint64_t         *il_Terms;
};

#endif

