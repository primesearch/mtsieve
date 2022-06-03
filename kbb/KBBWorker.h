/* KBBWorker.h -- (C) Mark Rodenkirch, July 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _KBBWorker_H
#define _KBBWorker_H

#include "KBBApp.h"
#include "../core/Worker.h"

using namespace std;

class KBBWorker : public Worker
{
public:
   KBBWorker(uint32_t myId, App *theApp);

   ~KBBWorker(void) {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   KBBApp           *ip_KBBApp;

   uint64_t          il_NextBaseBuild;
   uint32_t          ii_BaseCount;
   uint32_t         *ip_Bases;
   
   uint64_t          il_K;
};

#endif

