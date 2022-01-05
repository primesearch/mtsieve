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
   uint64_t  ps[4], maxPrime = ip_App->GetMaxPrime();
   vector<uint64_t>::iterator it = iv_Primes.begin();
   uint32_t *terms = ip_Terms->termList;
   uint32_t  termCount = ip_Terms->termCount;
   uint64_t  invmod2[4];
   uint64_t  invmod3[4];
   uint64_t  invmod4[4];
   bool      factorFound = false;
   
   uint64_t  m = ((uint64_t) terms[0] * 999999) + 1000000;
   uint32_t  exp = 6*terms[0] - 599989;

   while (it != iv_Primes.end())
   {
      ps[0] = *it;
      it++;
      
      ps[1] = *it;
      it++;
      
      ps[2] = *it;
      it++;
      
      ps[3] = *it;
      it++;

      invmod2[0] = InvMod64(12321, ps[0]);
      invmod2[1] = InvMod64(12321, ps[1]);
      invmod2[2] = InvMod64(12321, ps[2]);
      invmod2[3] = InvMod64(12321, ps[3]);

      invmod3[0] = InvMod64(1234321, ps[0]);
      invmod3[1] = InvMod64(1234321, ps[1]);
      invmod3[2] = InvMod64(1234321, ps[2]);
      invmod3[3] = InvMod64(1234321, ps[3]);

      invmod4[0] = InvMod64(123454321, ps[0]);
      invmod4[1] = InvMod64(123454321, ps[1]);
      invmod4[2] = InvMod64(123454321, ps[2]);
      invmod4[3] = InvMod64(123454321, ps[3]);
      
      MpArithVec mp(ps);

      MpResVec resTenE1 = mp.nToRes(10);
      
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
      
      MpResVec resTemp = mp.pow(resTenE1, 179);
      MpResVec resC = mp.mul(resTemp, tempMul1);
      resC = mp.sub(resC, tempSub1);

      //     t = powmod(10,2699,f);
      //     C = mulmod(C, t, f);
      //     submod(C,1109890222ULL);
      //     divmod(C, 12321ULL);

      resTemp = mp.pow(resTenE1, 2699);
      resC = mp.mul(resC, resTemp);
      resC = mp.sub(resC, tempSub2);
      resC = mp.mul(resC, resInvmod2);

      //     t = powmod(10,35999,f);
      //     C = mulmod(C, t, f);
      //     C = mulmod(C, 123454321ULL, f);
      //     submod(C,1110988902222ULL);
      //     divmod(C, 1234321ULL);
      
      resTemp = mp.pow(resTenE1, 35999);
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
      
      resTemp = mp.pow(resTenE1, 449999);
      resC = mp.mul(resC, resTemp);
      resC = mp.sub(resC, tempSub4);
      resC = mp.mul(resC, resInvmod4);
      
      //      t = powmod(10,6*v[0]-599989,f);
      //      C = mulmod(C, t, f);
      //      M = (v[0]*999999ULL+1000000)%f;
      //      r = 1000000000000000000ULL%f;
      
      resTemp = mp.pow(resTenE1, exp);
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
   
   if (factorFound)
   {
      xfree(ip_Terms);
      
      ip_Terms = ip_SmarandacheApp->GetTerms();
   }
}

void  SmarandacheWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("SmarandacheWorker::TestMiniPrimeChunk not implemented");
}
