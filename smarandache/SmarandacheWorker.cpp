/* SmarandacheWorker.cpp -- (C) Mark Rodenkirch, January 2022

   This program is free software;you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation;either version 2 of the License, or
   (at your option) any later version.
*/

#include <math.h>
#include "SmarandacheWorker.h"
#include "../core/MpArithVector.h"

extern "C" int mfsieve(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);
extern "C" int Smarandache(uint32_t start, uint32_t mf, uint32_t minmax, uint64_t *P);

SmarandacheWorker::SmarandacheWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_SmarandacheApp = (SmarandacheApp *) theApp;
   
   ii_MinN = ip_SmarandacheApp->GetMinN();
   ii_MaxN = ip_SmarandacheApp->GetMaxN();

   ip_Terms = ip_SmarandacheApp->GetTerms();
   
   ib_Initialized = true;
}

void  SmarandacheWorker::CleanUp(void)
{
   xfree(ip_Terms);
}

void  SmarandacheWorker::TestMegaPrimeChunk(void)
{
   uint32_t *terms = ip_Terms->termList;
   bool      factorFound;
   
   if (terms[0] < 1000000)
      factorFound = TestSixDigitN();
   else
      factorFound = TestSevenDigitN();
   
   if (factorFound)
   {
      xfree(ip_Terms);
      
      ip_Terms = ip_SmarandacheApp->GetTerms();
   }
}

bool  SmarandacheWorker::TestSixDigitN(void)
{
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   uint32_t *terms = ip_Terms->termList;
   uint32_t  termCount = ip_Terms->termCount;
   uint64_t  invmod2[4];
   uint64_t  invmod3[4];
   uint64_t  invmod4[4];
   bool      factorFound = false;
   
   uint64_t  m = ((uint64_t) terms[0] * 999999) + 1000000;
   uint32_t  exp = 6*terms[0] - 599989;
   
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      if (ps[0] < 12321)
      {
         invmod2[0] = InvMod64(12321 % ps[0], ps[0]);
         invmod2[1] = InvMod64(12321 % ps[1], ps[1]);
         invmod2[2] = InvMod64(12321 % ps[2], ps[2]);
         invmod2[3] = InvMod64(12321 % ps[3], ps[3]);
      }
      else
      {
         invmod2[0] = InvMod64(12321, ps[0]);
         invmod2[1] = InvMod64(12321, ps[1]);
         invmod2[2] = InvMod64(12321, ps[2]);
         invmod2[3] = InvMod64(12321, ps[3]);
      }

      if (ps[0] < 1234321)
      {
         invmod3[0] = InvMod64(1234321 % ps[0], ps[0]);
         invmod3[1] = InvMod64(1234321 % ps[1], ps[1]);
         invmod3[2] = InvMod64(1234321 % ps[2], ps[2]);
         invmod3[3] = InvMod64(1234321 % ps[3], ps[3]);
      }
      else
      {
         invmod3[0] = InvMod64(1234321, ps[0]);
         invmod3[1] = InvMod64(1234321, ps[1]);
         invmod3[2] = InvMod64(1234321, ps[2]);
         invmod3[3] = InvMod64(1234321, ps[3]);
      }
      
      if (ps[0] < 123454321)
      {
         invmod4[0] = InvMod64(123454321 % ps[0], ps[0]);
         invmod4[1] = InvMod64(123454321 % ps[1], ps[1]);
         invmod4[2] = InvMod64(123454321 % ps[2], ps[2]);
         invmod4[3] = InvMod64(123454321 % ps[3], ps[3]);
      }
      else
      {
         invmod4[0] = InvMod64(123454321, ps[0]);
         invmod4[1] = InvMod64(123454321, ps[1]);
         invmod4[2] = InvMod64(123454321, ps[2]);
         invmod4[3] = InvMod64(123454321, ps[3]);
      }
      
      
      MpArithVec mp(ps);

      MpResVec res10E1 = mp.nToRes(10);
      
      MpResVec tempMul1 = mp.nToRes(15208068915062105958ULL);
      MpResVec tempSub1 = mp.nToRes(11211123422ULL);
      
      MpResVec tempSub2 = mp.nToRes(1109890222);
      MpResVec resInvmod2 = mp.nToRes(invmod2);
      
      MpResVec tempMul3 = mp.nToRes(123454321);
      MpResVec tempSub3 = mp.nToRes(1110988902222ULL);
      MpResVec resInvmod3 = mp.nToRes(invmod3);
      
      MpResVec tempMul4 = mp.nToRes(12345654321ULL);
      MpResVec tempSub4 = mp.nToRes(1111098889022222ULL);
      MpResVec resInvmod4 = mp.nToRes(invmod4);
      
      //     t = powmod(10,179,f);
      //     C = mulmod(t,15208068915062105958ULL%f,f);
      //     submod(C,11211123422ULL);
      
      MpResVec resTemp = mp.pow(res10E1, 179);
      MpResVec resC = mp.mul(resTemp, tempMul1);
      resC = mp.sub(resC, tempSub1);

      //     t = powmod(10,2699,f);
      //     C = mulmod(C, t, f);
      //     submod(C,1109890222ULL);
      //     divmod(C, 12321ULL);

      resTemp = mp.pow(res10E1, 2699);
      resC = mp.mul(resC, resTemp);
      resC = mp.sub(resC, tempSub2);
      resC = mp.mul(resC, resInvmod2);

      //     t = powmod(10,35999,f);
      //     C = mulmod(C, t, f);
      //     C = mulmod(C, 123454321ULL, f);
      //     submod(C,1110988902222ULL);
      //     divmod(C, 1234321ULL);
      
      resTemp = mp.pow(res10E1, 35999);
      resC = mp.mul(resC, resTemp);
      resC = mp.mul(resC, tempMul3);
      resC = mp.sub(resC, tempSub3);
      resC = mp.mul(resC, resInvmod3);
      resC = mp.mul(resC, tempMul4);
      
      //      C = mulmod(C,12345654321ULL % f, f);
      //      t = powmod(10,449999,f);
      //      C = mulmod(C, t, f);
      //      submod(C,1111098889022222ULL);
      //      divmod(C, 123454321ULL);
      
      resTemp = mp.pow(res10E1, 449999);
      resC = mp.mul(resC, resTemp);
      resC = mp.sub(resC, tempSub4);
      resC = mp.mul(resC, resInvmod4);
      
      //      t = powmod(10,6*v[0]-599989,f);
      //      C = mulmod(C, t, f);
      //      M = (v[0]*999999ULL+1000000)%f;
      //      r = 1000000000000000000ULL%f;
      
      resTemp = mp.pow(res10E1, exp);
      resC = mp.mul(resC, resTemp);
      
      MpResVec resM = mp.nToRes(m);
      
      MpResVec res9s = mp.nToRes(999999);
      MpResVec resR = mp.nToRes(1000000000000000000ULL);
      MpResVec resT[6];

      //      T[1] = mulmod(r, r, f);
      //      T[2] = mulmod(T[1], T[1], f);
      //      T[3] = mulmod(T[1], T[2], f);
      //      T[4] = mulmod(T[2], T[2], f);
      //      T[5] = mulmod(T[2], T[3], f);

      resT[1] = mp.mul(resR, resR);
      resT[2] = mp.mul(resT[1], resT[1]);
      resT[3] = mp.mul(resT[1], resT[2]);
      resT[4] = mp.mul(resT[2], resT[2]);
      resT[5] = mp.mul(resT[2], resT[3]);
      
      // This is only for the first term
      if (MpArithVec::at_least_one_is_equal(resC, resM))
      {
         factorFound = true;
         
         for (uint32_t k=0;k<VECTOR_SIZE;++k)
            if (resC[k] == resM[k])
               ip_SmarandacheApp->ReportFactor(ps[k], terms[0]);
      }
         
      // This is for the remaining terms
      for (uint32_t idx=1; idx<termCount; idx++)
      {         
         uint32_t dn = terms[idx] - terms[idx-1];
         
         //    if(dn<=30)
         //      C = mulmod(C, T[dn/6], f);
         //    else {
         //      t = powmod(r, dn/3, f);
         //      C = mulmod(C, t, f);
         //    }
         
         if (dn <= 30)
            resC = mp.mul(resC, resT[dn/6]);
         else
         {
            resTemp = mp.pow(resT[1], dn/6);
            resC = mp.mul(resC, resTemp);
         }
         
         // M = (M + 999999ULL*dn) % f;
         
         resTemp = mp.nToRes(dn);
         resTemp = mp.mul(resTemp, res9s);
         resM = mp.add(resM, resTemp);

         if (MpArithVec::at_least_one_is_equal(resC, resM))
         {
            factorFound = true;
            
            for (uint32_t k=0;k<VECTOR_SIZE;++k)
               if (resC[k] == resM[k])
                  ip_SmarandacheApp->ReportFactor(ps[k], terms[idx]);
         }
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] >= maxPrime)
         break;
   }
   
   return factorFound;
}

bool  SmarandacheWorker::TestSevenDigitN(void)
{
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   uint32_t *terms = ip_Terms->termList;
   uint32_t  termCount = ip_Terms->termCount;
   uint64_t  invmod2[4], invmod3[4], invmod4[4], invmod5[4], invmod6[4];
   bool      factorFound = false;
   uint64_t  two9sq = 99;
   uint64_t  three9sq = 999;
   uint64_t  four9sq = 9999;
   uint64_t  five9sq = 99999;
   uint64_t  six9sq = 999999;
   uint64_t  seven9sq = 9999999;

   uint64_t  m = ((uint64_t) terms[0] * 9999999) + 10000000;
   uint32_t  exp = 7*terms[0] - 6999993;

   two9sq *= two9sq;
   three9sq *= three9sq;
   four9sq *= four9sq;
   five9sq *= five9sq;
   six9sq *= six9sq;
   seven9sq *= seven9sq;
   
   MpResVec resTemp;
   MpResVec resT[7];
   MpResVec resInvmod;
      
   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx+=4)
   {
      ps[0] = il_PrimeList[pIdx+0];
      ps[1] = il_PrimeList[pIdx+1];
      ps[2] = il_PrimeList[pIdx+2];
      ps[3] = il_PrimeList[pIdx+3];
      
      invmod2[0] = InvMod64(two9sq, ps[0]);
      invmod2[1] = InvMod64(two9sq, ps[1]);
      invmod2[2] = InvMod64(two9sq, ps[2]);
      invmod2[3] = InvMod64(two9sq, ps[3]);
     
      if (ps[0] < three9sq)
      {
         invmod3[0] = InvMod64(three9sq%ps[0], ps[0]);
         invmod3[1] = InvMod64(three9sq%ps[1], ps[1]);
         invmod3[2] = InvMod64(three9sq%ps[2], ps[2]);
         invmod3[3] = InvMod64(three9sq%ps[3], ps[3]);
      }
      else
      {
         invmod3[0] = InvMod64(three9sq, ps[0]);
         invmod3[1] = InvMod64(three9sq, ps[1]);
         invmod3[2] = InvMod64(three9sq, ps[2]);
         invmod3[3] = InvMod64(three9sq, ps[3]);
      }
      
      if (ps[0] < four9sq)
      {
         invmod4[0] = InvMod64(four9sq%ps[0], ps[0]);
         invmod4[1] = InvMod64(four9sq%ps[1], ps[1]);
         invmod4[2] = InvMod64(four9sq%ps[2], ps[2]);
         invmod4[3] = InvMod64(four9sq%ps[3], ps[3]);
      }
      else
      {
         invmod4[0] = InvMod64(four9sq, ps[0]);
         invmod4[1] = InvMod64(four9sq, ps[1]);
         invmod4[2] = InvMod64(four9sq, ps[2]);
         invmod4[3] = InvMod64(four9sq, ps[3]);
      }
      
      if (ps[0] < five9sq)
      {
         invmod5[0] = InvMod64(five9sq%ps[0], ps[0]);
         invmod5[1] = InvMod64(five9sq%ps[1], ps[1]);
         invmod5[2] = InvMod64(five9sq%ps[2], ps[2]);
         invmod5[3] = InvMod64(five9sq%ps[3], ps[3]);
      }
      else
      {
         invmod5[0] = InvMod64(five9sq, ps[0]);
         invmod5[1] = InvMod64(five9sq, ps[1]);
         invmod5[2] = InvMod64(five9sq, ps[2]);
         invmod5[3] = InvMod64(five9sq, ps[3]);
      }
      
      if (ps[0] < six9sq)
      {
         invmod6[0] = InvMod64(six9sq%ps[0], ps[0]);
         invmod6[1] = InvMod64(six9sq%ps[1], ps[1]);
         invmod6[2] = InvMod64(six9sq%ps[2], ps[2]);
         invmod6[3] = InvMod64(six9sq%ps[3], ps[3]);
      }
      else
      {
         invmod6[0] = InvMod64(six9sq, ps[0]);
         invmod6[1] = InvMod64(six9sq, ps[1]);
         invmod6[2] = InvMod64(six9sq, ps[2]);
         invmod6[3] = InvMod64(six9sq, ps[3]);
      }
      
      MpArithVec mp(ps);

      MpResVec res10e1 = mp.nToRes(10);
      MpResVec res10en = mp.nToRes(10);
      
      // calculate a1x
      //    ax = 123456789%f;
      MpResVec resAx = mp.nToRes(123456789);
      
      // calculate a2x
      //   ax = mulmod(ax, 99*99, f);
      //   ax = (ax + 991)%f;
      //    t = powmod(10,180,f);
      //   ax = mulmod(ax, t, f);
      //        submod(ax, 100);
      //        divmod(ax, 99*99);
      //        submod(ax, 1);
            
      resTemp = mp.nToRes(two9sq);
      resAx = mp.mul(resAx, resTemp);
      resTemp = mp.nToRes(991);
      resAx = mp.add(resAx, resTemp);
      resTemp = mp.pow(res10e1, 180);
      resAx = mp.mul(resAx, resTemp);
      res10en = mp.mul(res10e1, res10e1);
      resAx = mp.sub(resAx, res10en);
      resInvmod = mp.nToRes(invmod2);
      resAx = mp.mul(resAx, resInvmod);
      resAx = mp.sub(resAx, mp.one());

      // calculate a3x
      //   ax = mulmod(ax, 999*999, f);
      //   ax = (ax + 99901)%f;
      //    t = powmod(10,2700,f);
      //   ax = mulmod(ax, t, f);
      //        submod(ax, 1000);
      //        divmod(ax, 999*999);
      //        submod(ax, 1);
                  
      resTemp = mp.nToRes(three9sq);
      resAx = mp.mul(resAx, resTemp);
      resTemp = mp.nToRes(99901);
      resAx = mp.add(resAx, resTemp);
      resTemp = mp.pow(res10e1, 2700);
      resAx = mp.mul(resAx, resTemp);
      res10en = mp.mul(res10en, res10e1);
      resAx = mp.sub(resAx, res10en);
      resInvmod = mp.nToRes(invmod3);
      resAx = mp.mul(resAx, resInvmod);
      resAx = mp.sub(resAx, mp.one());
     
       // calculate a4x
      //   ax = mulmod(ax, 9999*9999, f);
      //   ax = (ax + 9999001)%f;
      //    t = powmod(10,36000,f);
      //   ax = mulmod(ax, t, f);
      //        submod(ax, 10000);
      //        divmod(ax, 9999*9999);
      //        submod(ax, 1);

      resTemp = mp.nToRes(four9sq);
      resAx = mp.mul(resAx, resTemp);
      resTemp = mp.nToRes(9999001);
      resAx = mp.add(resAx, resTemp);
      resTemp = mp.pow(res10e1, 36000);
      resAx = mp.mul(resAx, resTemp);
      res10en = mp.mul(res10en, res10e1);
      resAx = mp.sub(resAx, res10en);
      resInvmod = mp.nToRes(invmod4);
      resAx = mp.mul(resAx, resInvmod);
      resAx = mp.sub(resAx, mp.one());

      // calculate a5x
      //   ax = mulmod(ax, 99999ULL*99999ULL, f);
      //   ax = (ax + 999990001ULL)%f;
      //    t = powmod(10,450000,f);
      //   ax = mulmod(ax, t, f);
      //        submod(ax, 100000);
      //        divmod(ax, 99999ULL*99999ULL);
      //        submod(ax, 1);

      resTemp = mp.nToRes(five9sq);
      resAx = mp.mul(resAx, resTemp);
      resTemp = mp.nToRes(999990001ULL);
      resAx = mp.add(resAx, resTemp);
      resTemp = mp.pow(res10e1, 450000);
      resAx = mp.mul(resAx, resTemp);
      res10en = mp.mul(res10en, res10e1);
      resAx = mp.sub(resAx, res10en);
      resInvmod = mp.nToRes(invmod5);
      resAx = mp.mul(resAx, resInvmod);
      resAx = mp.sub(resAx, mp.one());

      // calculate a6x
      //   ax = mulmod(ax, 999999ULL*999999ULL, f);
      //   ax = (ax + 99999900001ULL)%f;
      //    t = powmod(10,5400000,f);
      //   ax = mulmod(ax, t, f);
      //        submod(ax, 1000000);
      //        divmod(ax, 999999ULL*999999ULL);
      //        submod(ax, 1);
         
      resTemp = mp.nToRes(six9sq);
      resAx = mp.mul(resAx, resTemp);
      resTemp = mp.nToRes(99999900001ULL);
      resAx = mp.add(resAx, resTemp);
      resTemp = mp.pow(res10e1, 5400000);
      resAx = mp.mul(resAx, resTemp);
      res10en = mp.mul(res10en, res10e1);
      resAx = mp.sub(resAx, res10en);
      resInvmod = mp.nToRes(invmod6);
      resAx = mp.mul(resAx, resInvmod);
      resAx = mp.sub(resAx, mp.one());
      
      // calculate a7(n)
      //   ax = mulmod(ax, 9999999 ULL*9999999ULL, f);
      //   ax = (ax + 9999999000001ULL)%f;
      //    t = powmod(10, 7*v[0]-6999993, f);
      //   ax = mulmod(ax, t, f);
      
      resTemp = mp.nToRes(seven9sq);
      resAx = mp.mul(resAx, resTemp);
      resTemp = mp.nToRes(9999999000001ULL);
      resAx = mp.add(resAx, resTemp);
      resTemp = mp.pow(res10e1, exp);
      resAx = mp.mul(resAx, resTemp);

      //    M = (9999999ULL*v[0]+10000000)%f;
      MpResVec resM = mp.nToRes(m);
      
      MpResVec res9s = mp.nToRes(9999999);
      MpResVec resTemp = mp.nToRes(10000000000ULL);
      MpResVec resR = mp.nToRes(100000000000ULL);

      //   r = mulmod(10000000000ULL, 100000000000ULL, f); // r = 10^21%f
      //   T[1] = mulmod(r, r, f);
      //   T[2] = mulmod(T[1], T[1], f);
      //   T[3] = mulmod(T[1], T[2], f);
      //   T[4] = mulmod(T[2], T[2], f);
      //   T[5] = mulmod(T[2], T[3], f);
      //   T[6] = mulmod(T[3], T[3], f);

      resR = mp.mul(resR, resTemp);
      resT[1] = mp.mul(resR, resR);
      resT[2] = mp.mul(resT[1], resT[1]);
      resT[3] = mp.mul(resT[1], resT[2]);
      resT[4] = mp.mul(resT[2], resT[2]);
      resT[5] = mp.mul(resT[2], resT[3]);
      resT[6] = mp.mul(resT[3], resT[3]);

      // This is only for the first term
      if (MpArithVec::at_least_one_is_equal(resAx, resM))
      {
         factorFound = true;
         
         for (uint32_t k=0;k<VECTOR_SIZE;++k)
            if (resAx[k] == resM[k])
               ip_SmarandacheApp->ReportFactor(ps[k], terms[0]);
      }
         
      // This is for the remaining terms
      for (uint32_t idx=1; idx<termCount; idx++)
      {         
         uint32_t dn = terms[idx] - terms[idx-1];
         
         //   if (dn<=36)
         //      ax = mulmod(ax, T[dn/6], f);
         //   else {
         //     t = powmod(T[1], dn/6, f);
         //     ax = mulmod(ax, t, f);
         //   }
         
         if (dn <= 36)
            resAx = mp.mul(resAx, resT[dn/6]);
         else
         {
            resTemp = mp.pow(resT[1], dn/6);
            resAx = mp.mul(resAx, resTemp);
         }
         
         // M = (M + 9999999ULL*dn) % f;
         
         resTemp = mp.nToRes(dn);
         resTemp = mp.mul(resTemp, res9s);
         resM = mp.add(resM, resTemp);

         if (MpArithVec::at_least_one_is_equal(resAx, resM))
         {
            factorFound = true;
            
            for (uint32_t k=0;k<VECTOR_SIZE;++k)
               if (resAx[k] == resM[k])
                  ip_SmarandacheApp->ReportFactor(ps[k], terms[idx]);
         }
      }
      
      SetLargestPrimeTested(ps[3], 4);
      
      if (ps[3] >= maxPrime)
         break;
   }
   
   return factorFound;
}

void  SmarandacheWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheWorker::TestMiniPrimeChunk not implemented");
}
