#ifdef __cplusplus
#ifndef _CISONESINGLE_KERNEL_CL
#define _CISONESINGLE_KERNEL_CL
const char *cisonesingle_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"#define SP_MIXED     0\n" \
"#define SP_EVEN      1\n" \
"#define SP_ODD       2\n" \
"#define SP_NO_PARITY 999\n" \
"#define N_TERM(q, i, j)       ((SIEVE_LOW + (j) + (i)*bSteps)*BESTQ + q)\n" \
"#define CSS_INDEX(x, y, z)    (((((x) * (PRL_COUNT + 1)) + (y)) * (POWER_RESIDUE_LCM + 1)) + (z))\n" \
"#define HASH_NOT_FOUND        UINT_MAX\n" \
"#define HASH_MASK1            (1<<15)\n" \
"#define HASH_MASK2            (HASH_MASK1-1)\n" \
"ushort getParity(ulong thePrime);\n" \
"ulong  invmod(ulong a, ulong p);\n" \
"short  legendre(long a, ulong p);\n" \
"uint  setupDiscreteLog(ulong thePrime, ulong _q, ulong _one,\n" \
"ulong resBase, ulong resInvBase, ulong resNegCK, ushort parity,\n" \
"__global const short  *divisorShifts,\n" \
"__global const uint   *prlIndices);\n" \
"ulong buildLookupsAndClimbLadder(ulong thePrime, ulong _q, ulong _one,\n" \
"ulong resBase, ulong resNegCK, ulong *resBD,\n" \
"uint qIndex, uint ladderIndex,\n" \
"__global const ushort *qs,\n" \
"__global const ushort *ladders);\n" \
"ulong babySteps(ulong thePrime, ulong _q, ulong _one, ulong resInvBase, uint bSteps,\n" \
"ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"void  hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"uint  hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"void collectFactor(ulong   p,\n" \
"uint    n,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong2 *factors);\n" \
"ulong mmmInvert(ulong p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);\n" \
"__kernel void cisonesingle_kernel(__global const ulong  *primes,\n" \
"__global const uint   *babyStepsArray,\n" \
"__global const uint   *giantStepsArray,\n" \
"__global const short  *divisorShifts,\n" \
"__global const uint   *prlIndices,\n" \
"__global const uint   *qIndices,\n" \
"__global const ushort *qs,\n" \
"__global const uint   *ladderIndices,\n" \
"__global const ushort *ladders,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong2 *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"ulong  thePrime = primes[gid];\n" \
"ulong  negCK;\n" \
"ushort parity;\n" \
"parity = getParity(thePrime);\n" \
"if (parity == SP_NO_PARITY)\n" \
"return;\n" \
"if (thePrime < SEQ_K)\n" \
"negCK = SEQ_K - (SEQ_K / thePrime);\n" \
"else\n" \
"negCK = SEQ_K;\n" \
"if (SEQ_C > 0)\n" \
"negCK = thePrime - negCK;\n" \
"ulong   _q = mmmInvert(thePrime);\n" \
"ulong   _one = mmmOne(thePrime);\n" \
"ulong   _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong   invBase = invmod(BASE, thePrime);\n" \
"ulong   resBase = mmmNToRes(BASE, thePrime, _q, _r2);\n" \
"ulong   resInvBase = mmmNToRes(invBase, thePrime, _q, _r2);\n" \
"ulong   resNegCK = mmmNToRes(negCK, thePrime, _q, _r2);\n" \
"ulong   resBexpQ;\n" \
"ushort  h_table[HASH_SIZE];\n" \
"ushort  h_olist[HASH_ELEMENTS];\n" \
"ulong   h_BJ64[HASH_ELEMENTS+1];\n" \
"ulong   resBD[SUBSEQUENCE_COUNT];\n" \
"uint    i, j, k, ssCount;\n" \
"uint    idx, qIdx, ladderIdx;\n" \
"ushort cssIndex = setupDiscreteLog(thePrime, _q, _one, resBase, resInvBase, resNegCK, parity, divisorShifts, prlIndices);\n" \
"qIdx = qIndices[cssIndex];\n" \
"if (qIdx == 0)\n" \
"return;\n" \
"ladderIdx = ladderIndices[cssIndex];\n" \
"resBexpQ = buildLookupsAndClimbLadder(thePrime, _q, _one, resBase, resNegCK, resBD, qIdx, ladderIdx, qs, ladders);\n" \
"ssCount = qs[qIdx];\n" \
"for (idx=0; idx<HASH_SIZE; idx++)\n" \
"h_table[idx] = 0;\n" \
"for (idx=0; idx<HASH_ELEMENTS; idx++)\n" \
"h_BJ64[idx] = 0;\n" \
"h_BJ64[HASH_ELEMENTS] = ULONG_MAX;\n" \
"uint bSteps = babyStepsArray[ssCount - 1];\n" \
"uint gSteps = giantStepsArray[ssCount - 1];\n" \
"ushort orderOfB = babySteps(thePrime, _q, _one, resInvBase, bSteps, h_table, h_olist, h_BJ64);\n" \
"if (orderOfB > 0)\n" \
"{\n" \
"for (k=0; k<ssCount; k++)\n" \
"{\n" \
"j = hashLookup(resBD[k], h_table, h_olist, h_BJ64);\n" \
"while (j < bSteps * gSteps)\n" \
"{\n" \
"collectFactor(thePrime, N_TERM(qs[qIdx+1+k], 0, j), factorCount, factors);\n" \
"j += orderOfB;\n" \
"}\n" \
"}\n" \
"}\n" \
"else\n" \
"{\n" \
"for (k=0; k<ssCount; k++)\n" \
"{\n" \
"j = hashLookup(resBD[k], h_table, h_olist, h_BJ64);\n" \
"if (j != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, N_TERM(qs[qIdx+1+k], 0, j), factorCount, factors);\n" \
"}\n" \
"if (gSteps > 1)\n" \
"{\n" \
"ulong resBQM = mmmPowmod(resBexpQ, bSteps, thePrime, _q, _one);\n" \
"for (i=1; i<gSteps; i++)\n" \
"{\n" \
"for (k=0; k<ssCount; k++)\n" \
"{\n" \
"resBD[k] = mmmMulmod(resBD[k], resBQM, thePrime, _q);\n" \
"j = hashLookup(resBD[k], h_table, h_olist, h_BJ64);\n" \
"if (j != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, N_TERM(qs[qIdx+1+k], i, j), factorCount, factors);\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"}\n" \
"ulong mmmInvert(ulong p)\n" \
"{\n" \
"ulong p_inv = 1;\n" \
"ulong prev = 0;\n" \
"while (p_inv != prev)\n" \
"{\n" \
"prev = p_inv;\n" \
"p_inv *= (2 - p * p_inv);\n" \
"}\n" \
"return p_inv;\n" \
"}\n" \
"inline ulong mmmOne(ulong _p)\n" \
"{\n" \
"return ((-_p) % _p);\n" \
"}\n" \
"inline ulong mmmR2(ulong _p, ulong _q, ulong _one)\n" \
"{\n" \
"ulong t = mmmAdd(_one, _one, _p);\n" \
"t = mmmAdd(t, t, _p);\n" \
"for (size_t i=0; i<5; i++)\n" \
"t = mmmMulmod(t, t, _p, _q);\n" \
"return t;\n" \
"}\n" \
"inline ulong mmmAdd(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a >= _p - b) ? _p : 0;\n" \
"return a + b - c;\n" \
"}\n" \
"inline ulong mmmSub(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a < b) ? _p : 0;\n" \
"return a - b + c;\n" \
"}\n" \
"inline ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
"}\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q)\n" \
"{\n" \
"ulong lo = a * b;\n" \
"ulong hi = mul_hi(a, b);\n" \
"ulong m = lo * _q;\n" \
"ulong hi2 = mul_hi(m, _p);\n" \
"long r = (long) hi - (long) hi2;\n" \
"if (r < 0)\n" \
"return (ulong) (r + _p);\n" \
"return (ulong) r;\n" \
"}\n" \
"ulong   mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one)\n" \
"{\n" \
"ulong x = resbase;\n" \
"ulong y = _one;\n" \
"while (true)\n" \
"{\n" \
"if (exp & 1)\n" \
"y = mmmMulmod(x, y, _p, _q);\n" \
"exp >>= 1;\n" \
"if (!exp)\n" \
"return y;\n" \
"x = mmmMulmod(x, x, _p, _q);\n" \
"}\n" \
"return 0;\n" \
"}\n" \
"ushort  getParity(ulong thePrime)\n" \
"{\n" \
"if (SEQ_PARITY == SP_MIXED)\n" \
"{\n" \
"short qr_m1, qr_p1;\n" \
"{\n" \
"short sym = legendre(KC_CORE, thePrime);\n" \
"qr_m1 = (sym == 1);\n" \
"qr_p1 = (sym == legendre(BASE, thePrime));\n" \
"}\n" \
"if (qr_m1)\n" \
"return (qr_p1 ? SP_MIXED : SP_EVEN);\n" \
"return (qr_p1 ? SP_ODD : SP_NO_PARITY);\n" \
"}\n" \
"short qr;\n" \
"{\n" \
"short sym = legendre(KC_CORE, thePrime);\n" \
"if (SEQ_PARITY == SP_EVEN)\n" \
"qr = (sym == 1);\n" \
"else\n" \
"qr = (sym == legendre(BASE, thePrime));\n" \
"}\n" \
"if (qr)\n" \
"return SEQ_PARITY;\n" \
"return SP_NO_PARITY;\n" \
"}\n" \
"uint setupDiscreteLog(ulong thePrime, ulong _q, ulong _one,\n" \
"ulong resBase, ulong resInvBase, ulong resNegCK, ushort parity,\n" \
"__global const short  *divisorShifts,\n" \
"__global const uint   *prlIndices)\n" \
"{\n" \
"ulong pShift;\n" \
"uint  idx;\n" \
"uint  h, r;\n" \
"short shift;\n" \
"ulong resX[POWER_RESIDUE_LCM + 1];\n" \
"idx = (thePrime/2) % (POWER_RESIDUE_LCM/2);\n" \
"shift = divisorShifts[idx];\n" \
"if (shift == 0)\n" \
"{\n" \
"r = prlIndices[1];\n" \
"h = 0;\n" \
"return CSS_INDEX(parity, r, h);\n" \
"}\n" \
"if (shift > 0)\n" \
"{\n" \
"pShift = thePrime / shift;\n" \
"}\n" \
"else\n" \
"{\n" \
"pShift = thePrime >> (-shift);\n" \
"shift = 1 << (-shift);\n" \
"}\n" \
"resX[0] = _one;\n" \
"resX[1] = mmmPowmod(resInvBase, pShift, thePrime, _q, _one);\n" \
"for (r=1; resX[r] != resX[0]; r++)\n" \
"resX[r+1] = mmmMulmod(resX[r], resX[1], thePrime, _q);\n" \
"if (shift % r != 0)\n" \
"return 0;\n" \
"resX[r] = mmmPowmod(resNegCK, pShift, thePrime, _q, _one);\n" \
"for (h=0; resX[r] != resX[h]; h++)\n" \
";\n" \
"if (h == r)\n" \
"return 0;\n" \
"r = prlIndices[r];\n" \
"return CSS_INDEX(parity, r, h);\n" \
"}\n" \
"ulong buildLookupsAndClimbLadder(ulong thePrime, ulong _q,  ulong _one,\n" \
"ulong resBase, ulong resNegCK, ulong *resBD,\n" \
"uint qIndex, uint ladderIndex,\n" \
"__global const ushort *qs,\n" \
"__global const ushort *ladders)\n" \
"{\n" \
"uint   i, j, idx, lLen, qLen;\n" \
"ulong  resX[POWER_RESIDUE_LCM + 1];\n" \
"lLen = ladders[ladderIndex];\n" \
"resX[0] = _one;\n" \
"resX[1] = resBase;\n" \
"resX[2] = mmmMulmod(resX[1], resX[1], thePrime, _q);\n" \
"i = 2;\n" \
"for (j=0; j<lLen; j++)\n" \
"{\n" \
"idx = ladders[ladderIndex+j+1];\n" \
"resX[i+idx] = mmmMulmod(resX[i], resX[idx], thePrime, _q);\n" \
"i += idx;\n" \
"}\n" \
"qLen = qs[qIndex];\n" \
"for (j=0; j<qLen; j++)\n" \
"{\n" \
"idx = qs[qIndex+j+1];\n" \
"resBD[j] = mmmMulmod(resX[idx], resNegCK, thePrime, _q);\n" \
"}\n" \
"return resX[BESTQ];\n" \
"}\n" \
"ulong babySteps(ulong thePrime, ulong _q, ulong _one,\n" \
"ulong resInvBase, uint bSteps,\n" \
"ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"ulong    resInvBaseExpQ;\n" \
"ulong    resBJ;\n" \
"ulong    firstResBJ;\n" \
"uint     j;\n" \
"resInvBaseExpQ = mmmPowmod(resInvBase, BESTQ, thePrime, _q, _one);\n" \
"firstResBJ = resBJ = mmmPowmod(resInvBaseExpQ, SIEVE_LOW, thePrime, _q, _one);\n" \
"for (j=0; j<bSteps; j++)\n" \
"{\n" \
"hashInsert(resBJ, j, h_table, h_olist, h_BJ64);\n" \
"resBJ = mmmMulmod(resBJ, resInvBaseExpQ, thePrime, _q);\n" \
"if (resBJ == firstResBJ)\n" \
"return j + 1;\n" \
"}\n" \
"return 0;\n" \
"}\n" \
"ulong  invmod(ulong a, ulong p)\n" \
"{\n" \
"ulong ps1, ps2, q, r, t, dividend, divisor;\n" \
"uint parity;\n" \
"if (a < 3)\n" \
"return (a < 2) ? a : (p+1)/2;\n" \
"q = p / a;\n" \
"r = p % a;\n" \
"dividend = a;\n" \
"divisor = r;\n" \
"ps1 = q;\n" \
"ps2 = 1;\n" \
"parity = 0;\n" \
"while (divisor > 1)\n" \
"{\n" \
"r = dividend - divisor;\n" \
"t = r - divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t; t -= divisor;\n" \
"if (r >= divisor) {\n" \
"q += ps1; r = t;\n" \
"if (r >= divisor) {\n" \
"q = dividend / divisor;\n" \
"r = dividend % divisor;\n" \
"q *= ps1;\n" \
"} } } } } } } } }\n" \
"q += ps2;\n" \
"parity = ~parity;\n" \
"dividend = divisor;\n" \
"divisor = r;\n" \
"ps2 = ps1;\n" \
"ps1 = q;\n" \
"}\n" \
"return (parity) ? ps1 : p - ps1;\n" \
"}\n" \
"short  legendre(long a, ulong p)\n" \
"{\n" \
"ulong x, y, t;\n" \
"short sign;\n" \
"if (a < 0)\n" \
"{\n" \
"a = -a;\n" \
"sign = (p % 4 == 1) ? 1 : -1;\n" \
"}\n" \
"else\n" \
"sign = 1;\n" \
"for (y = a; y % 2 == 0; y /= 2)\n" \
"if (p % 8 == 3 || p % 8 == 5)\n" \
"sign = -sign;\n" \
"if (p % 4 == 3 && y % 4 == 3)\n" \
"sign = -sign;\n" \
"for (x = p % y; x > 0; x %= y)\n" \
"{\n" \
"for ( ; x % 2 == 0; x /= 2)\n" \
"if (y % 8 == 3 || y % 8 == 5)\n" \
"sign = -sign;\n" \
"t = x, x = y, y = t;\n" \
"if (x % 4 == 3 && y % 4 == 3)\n" \
"sign = -sign;\n" \
"}\n" \
"return sign;\n" \
"}\n" \
"void hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint slot;\n" \
"h_BJ64[j] = bj;\n" \
"slot = bj & (HASH_SIZE - 1);\n" \
"if (h_table[slot] == HASH_ELEMENTS)\n" \
"h_table[slot] = j;\n" \
"else\n" \
"{\n" \
"h_olist[j] = (h_table[slot] ^ HASH_MASK1);\n" \
"h_table[slot] = (j | HASH_MASK1);\n" \
"}\n" \
"}\n" \
"uint hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint   slot;\n" \
"ushort elt;\n" \
"slot = bj & (HASH_SIZE - 1);\n" \
"elt = h_table[slot];\n" \
"if (h_BJ64[elt & HASH_MASK2] == bj)\n" \
"return elt & HASH_MASK2;\n" \
"if ((elt & HASH_MASK1) == 0)\n" \
"return HASH_NOT_FOUND;\n" \
"elt &= HASH_MASK2;\n" \
"do\n" \
"{\n" \
"if (h_BJ64[h_olist[elt] & HASH_MASK2] == bj)\n" \
"return h_olist[elt] & HASH_MASK2;\n" \
"elt = h_olist[elt];\n" \
"} while ((elt & HASH_MASK1) == 0);\n" \
"return HASH_NOT_FOUND;\n" \
"}\n" \
"void collectFactor(ulong   p,\n" \
"uint    n,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong2 *factors)\n" \
"{\n" \
"int old = atomic_inc(factorCount);\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = n;\n" \
"factors[old].y = p;\n" \
"}\n" \
;
#endif
#endif
