/* AbstractWorker.h -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _AbstractWorker_H
#define _AbstractWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSubsequenceHelper.h"
#include "../core/Worker.h"
#include "../core/HashTable.h"

using namespace std;

class AbstractWorker : public Worker
{
public:
   AbstractWorker(uint32_t myId, App *theApp, AbstractSubsequenceHelper *appHelper);

   ~AbstractWorker(void) {};
   
   virtual void         SetSequences(uint64_t largestPrimeTested, uint32_t bestQ, seq_t *sequences, uint32_t sequenceCount, subseq_t *subsequences, uint32_t subsequenceCount) = 0;
   
protected:
   virtual void         InitializeWorker(void) = 0;
   
   SierpinskiRieselApp *ip_SierpinskiRieselApp;
   AbstractSubsequenceHelper   *ip_AppHelper;
   
   uint32_t             ii_Base;
   uint32_t             ii_MinN;
   uint32_t             ii_MaxN;
   
   // The sequences and subsequences are read only in the worker classes
   seq_t               *ip_Sequences;
   subseq_t            *ip_Subsequences;
   
   uint32_t             ii_SequenceCount;
   uint32_t             ii_SubsequenceCount;

   uint32_t             ii_BestQ;
   
   // 1/a (mod p)  if a > 0, 0 othewise. Assumes a < p and gcd(a,p)=1.
   // Thanks to the folks at mersenneforum.org.
   // See http://www.mersenneforum.org/showthread.php?p=58252.
   inline uint64_t  invmod64(uint64_t a, uint64_t p)
   {
      uint64_t ps1, ps2, q, r, t, dividend, divisor;
      uint32_t parity;

      if (a >= p)
         FatalError("invmod called with invalid parameters (%" PRIu64" > %" PRIu64")", a, p);

      if (a < 3)
         return (a < 2) ? a : (p+1)/2;

      q = p / a;
      r = p % a;

      if (r == 0)
         FatalError("invmod called with a % p = 0");

      dividend = a;
      divisor = r;
      ps1 = q;
      ps2 = 1;
      parity = 0;

      while (divisor > 1)
      {
         r = dividend - divisor;
         t = r - divisor;
         if (r >= divisor) {
            q += ps1; r = t; t -= divisor;
            if (r >= divisor) {
               q += ps1; r = t; t -= divisor;
               if (r >= divisor) {
                  q += ps1; r = t; t -= divisor;
                  if (r >= divisor) {
                     q += ps1; r = t; t -= divisor;
                     if (r >= divisor) {
                        q += ps1; r = t; t -= divisor;
                        if (r >= divisor) {
                           q += ps1; r = t; t -= divisor;
                           if (r >= divisor) {
                              q += ps1; r = t; t -= divisor;
                              if (r >= divisor) {
                                 q += ps1; r = t;
                                 if (r >= divisor) {
                                    q = dividend / divisor;
                                    r = dividend % divisor;
                                    q *= ps1;
                                 } } } } } } } } }
         q += ps2;
         parity = ~parity;
         dividend = divisor;
         divisor = r;
         ps2 = ps1;
         ps1 = q;
      }

      if (ps1 == 0)
         FatalError("invmod failed with ps1 = 0");
      
      if (ps1 == 0)
         FatalError("invmod failed with ps1 >= p");

      return (parity) ? ps1 : p - ps1;
   }

   // return a (mod p), assuming -p < a < p
   inline uint64_t  lmod64(int64_t a, uint64_t p)
   {
      uint64_t ua;
      if (a >= 0)
      {
         ua = (uint64_t) a;
         
         return (ua < p) ? ua : ua%p;
      }
      else
      {
         ua = (uint64_t) -a;
         
         return (ua < p) ? p-ua : p-(ua%p);
      }
   }

   // return a (mod p), assuming 0 <= a
   inline uint64_t  umod64(uint64_t a, uint64_t p)
   {
      return (a < p) ? a : a%p;
   }
   
private:

};

#endif

