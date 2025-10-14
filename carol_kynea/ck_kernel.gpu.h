#ifdef __cplusplus
#ifndef _CK_KERNEL_CL
#define _CK_KERNEL_CL
const char *ck_kernel= \
"#if !defined(USE_OPENCL) && !defined(USE_METAL)\n" \
"#define MEED_METALLIB\n" \
"#endif\n" \
"#ifdef USE_METAL\n" \
"#include <metal_stdlib>\n" \
"#include <metal_atomic>\n" \
"using namespace metal;\n" \
"#endif\n" \
"#ifdef MEED_METALLIB\n" \
"#define BASE                 2\n" \
"#define SIEVE_LOW           50\n" \
"#define SIEVE_RANGE         50\n" \
"#define BABY_STEPS          50\n" \
"#define GIANT_STEPS         50\n" \
"#define HASH_ELEMENTS      500\n" \
"#define HASH_SIZE          500\n" \
"#define MAX_FACTORS        500\n" \
"#endif\n" \
"#if defined(USE_OPENCL)\n" \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MUL_HI          mul_hi\n" \
"void collectFactor(ulong   p,\n" \
"ulong   n,\n" \
"int     c,\n" \
"volatile __global uint   *factorCount,\n" \
"__global long4  *factors);\n" \
"#else\n" \
"#define MUL_HI          mulhi\n" \
"void collectFactor(ulong       p,\n" \
"ulong       n,\n" \
"int         c,\n" \
"volatile device atomic_int *factorCount,\n" \
"device long4      *factors);\n" \
"#endif\n" \
"#define HASH_NOT_FOUND        UINT_MAX\n" \
"#define HASH_MASK1            (1<<15)\n" \
"#define HASH_MASK2            (HASH_MASK1-1)\n" \
"#define ROOT_COUNT            4\n" \
"ulong invmod(ulong a, ulong p);\n" \
"uint  isQuadraticResidue(ulong n, ulong p);\n" \
"ulong findRoot(ulong thePrime, ulong _q, ulong _r2, ulong _one);\n" \
"ulong babySteps(ulong thePrime, ulong _q, ulong _r2, ulong _one, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"void  hashInsert(ulong bj, uint j, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"uint  hashLookup(ulong bj, ushort *h_table, ushort *h_olist, ulong *h_BJ64);\n" \
"ulong mmmInvert(ulong p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);\n" \
"#if defined(USE_OPENCL)\n" \
"__kernel void ck_kernel(__global const ulong  *primes,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       long4  *factors)\n" \
"#else\n" \
"kernel void ck_kernel(device const ulong  *primes,\n" \
"volatile device atomic_int   *factorCount,\n" \
"device       long4  *factors\n" \
"uint    tpid [[thread_position_in_grid]])\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int    gid = get_global_id(0);\n" \
"#else\n" \
"int    gid = tpid;\n" \
"#endif\n" \
"ushort  h_table[HASH_SIZE];\n" \
"ushort  h_olist[HASH_ELEMENTS];\n" \
"ulong   h_BJ64[HASH_ELEMENTS+1];\n" \
"uint    jHash[ROOT_COUNT];\n" \
"ulong   resA[ROOT_COUNT];\n" \
"uint    j, idx;\n" \
"ulong   root1, root2;\n" \
"ulong   thePrime = primes[gid];\n" \
"if (BASE % thePrime == 0)\n" \
"return;\n" \
"if (isQuadraticResidue(2, thePrime) == 0)\n" \
"return;\n" \
"ulong inv_pb = invmod(BASE, thePrime);\n" \
"if (inv_pb == 0)\n" \
"return;\n" \
"for (idx=0; idx<HASH_SIZE; idx++)\n" \
"h_table[idx] = HASH_ELEMENTS;\n" \
"for (idx=0; idx<HASH_ELEMENTS; idx++)\n" \
"h_BJ64[idx] = 0;\n" \
"h_BJ64[HASH_ELEMENTS] = ULONG_MAX;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"root1 = findRoot(thePrime, _q, _r2, _one);\n" \
"root2 = mmmSub(mmmNToRes(thePrime, thePrime, _q, _r2), root1, thePrime);\n" \
"resA[0] = mmmSub(root1, _one, thePrime);\n" \
"resA[1] = mmmSub(root2, _one, thePrime);\n" \
"resA[2] = mmmAdd(root1, _one, thePrime);\n" \
"resA[3] = mmmAdd(root2, _one, thePrime);\n" \
"uint orderOfB = babySteps(thePrime, _q, _r2, _one, h_table, h_olist, h_BJ64);\n" \
"jHash[0] = hashLookup(resA[0], h_table, h_olist, h_BJ64);\n" \
"jHash[1] = hashLookup(resA[1], h_table, h_olist, h_BJ64);\n" \
"jHash[2] = hashLookup(resA[2], h_table, h_olist, h_BJ64);\n" \
"jHash[3] = hashLookup(resA[3], h_table, h_olist, h_BJ64);\n" \
"if (orderOfB > 0)\n" \
"{\n" \
"uint j;\n" \
"for (j = jHash[0]; j < SIEVE_RANGE; j += orderOfB)\n" \
"collectFactor(thePrime, SIEVE_LOW+j, +1, factorCount, factors);\n" \
"for (j = jHash[1]; j < SIEVE_RANGE; j += orderOfB)\n" \
"collectFactor(thePrime, SIEVE_LOW+j, +1, factorCount, factors);\n" \
"for (j = jHash[2]; j < SIEVE_RANGE; j += orderOfB)\n" \
"collectFactor(thePrime, SIEVE_LOW+j, -1, factorCount, factors);\n" \
"for (j = jHash[3]; j < SIEVE_RANGE; j += orderOfB)\n" \
"collectFactor(thePrime, SIEVE_LOW+j, -1, factorCount, factors);\n" \
"return;\n" \
"}\n" \
"if (jHash[0] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, SIEVE_LOW+jHash[0], +1, factorCount, factors);\n" \
"if (jHash[1] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, SIEVE_LOW+jHash[1], +1, factorCount, factors);\n" \
"if (jHash[2] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, SIEVE_LOW+jHash[2], -1, factorCount, factors);\n" \
"if (jHash[3] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, SIEVE_LOW+jHash[3], -1, factorCount, factors);\n" \
"ulong resB = mmmNToRes(inv_pb, thePrime, _q, _r2);\n" \
"resB = mmmPowmod(resB, BABY_STEPS, thePrime, _q, _one);\n" \
"uint   nBase = SIEVE_LOW;\n" \
"for (uint step=1; step<GIANT_STEPS; step++)\n" \
"{\n" \
"resA[0] = mmmMulmod(resA[0], resB, thePrime, _q);\n" \
"resA[1] = mmmMulmod(resA[1], resB, thePrime, _q);\n" \
"resA[2] = mmmMulmod(resA[2], resB, thePrime, _q);\n" \
"resA[3] = mmmMulmod(resA[3], resB, thePrime, _q);\n" \
"jHash[0] = hashLookup(resA[0], h_table, h_olist, h_BJ64);\n" \
"jHash[1] = hashLookup(resA[1], h_table, h_olist, h_BJ64);\n" \
"jHash[2] = hashLookup(resA[2], h_table, h_olist, h_BJ64);\n" \
"jHash[3] = hashLookup(resA[3], h_table, h_olist, h_BJ64);\n" \
"nBase += BABY_STEPS;\n" \
"if (jHash[0] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, nBase+jHash[0], +1, factorCount, factors);\n" \
"if (jHash[1] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, nBase+jHash[1], +1, factorCount, factors);\n" \
"if (jHash[2] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, nBase+jHash[2], -1, factorCount, factors);\n" \
"if (jHash[3] != HASH_NOT_FOUND)\n" \
"collectFactor(thePrime, nBase+jHash[3], -1, factorCount, factors);\n" \
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
"ulong hi = MUL_HI(a, b);\n" \
"ulong m = lo * _q;\n" \
"ulong hi2 = MUL_HI(m, _p);\n" \
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
"ulong mmmResToN(ulong res, ulong _p, ulong _q)\n" \
"{\n" \
"return mmmMulmod(res, 1, _p, _q);\n" \
"}\n" \
"uint  isQuadraticResidue(ulong n, ulong _p)\n" \
"{\n" \
"int      j = 1;\n" \
"ulong    swap;\n" \
"while (n > 1)\n" \
"{\n" \
"if (n & 1)\n" \
"{\n" \
"swap = n;\n" \
"n = _p;\n" \
"_p = swap;\n" \
"if (((n & 3) == 3) && ((_p & 3) == 3))\n" \
"j = -j;\n" \
"n = n % _p;\n" \
"}\n" \
"else\n" \
"{\n" \
"n >>= 1;\n" \
"if (((_p & 7) == 3) || ((_p & 7) == 5))\n" \
"j = -j;\n" \
"}\n" \
"}\n" \
"if (j == 1)\n" \
"return 1;\n" \
"return 0;\n" \
"}\n" \
"ulong findRoot(ulong thePrime, ulong _q, ulong _r2, ulong _one)\n" \
"{\n" \
"ulong	   i, s, t, d, m;\n" \
"ulong    res, resA, resD;\n" \
"ulong    res2 = mmmNToRes(2, thePrime, _q, _r2);\n" \
"ulong    resPM1 = mmmNToRes(thePrime-1, thePrime, _q, _r2);\n" \
"if ((thePrime & 7) == 7)\n" \
"return mmmPowmod(res2, (thePrime+1) >> 2, thePrime, _q, _one);\n" \
"t = thePrime - 1;\n" \
"s = 0;\n" \
"while (!(t & 1))\n" \
"{\n" \
"s++;\n" \
"t >>= 1;\n" \
"}\n" \
"resA =  mmmPowmod(res2, t, thePrime, _q, _one);\n" \
"for (d=3; d<thePrime; d++)\n" \
"if (!isQuadraticResidue(d, thePrime))\n" \
"break;\n" \
"resD = mmmNToRes(d, thePrime, _q, _r2);\n" \
"resD = mmmPowmod(resD, t, thePrime, _q, _one);\n" \
"m = 0;\n" \
"for (i=0; i<s; i++)\n" \
"{\n" \
"if (m == 0)\n" \
"res = _one;\n" \
"else\n" \
"res = mmmPowmod(resD, m, thePrime, _q, _one);\n" \
"res = mmmMulmod(res, resA, thePrime, _q);\n" \
"res = mmmPowmod(res, 1 << (s - 1 - i), thePrime, _q, _one);\n" \
"if (res == resPM1)\n" \
"m += (1 << i);\n" \
"}\n" \
"if (m == 0)\n" \
"res = _one;\n" \
"else\n" \
"res = mmmPowmod(resD, m >> 1, thePrime, _q, _one);\n" \
"resA = mmmPowmod(res2, (t+1) >> 1, thePrime, _q, _one);\n" \
"return mmmMulmod(res, resA, thePrime, _q);\n" \
"}\n" \
"ulong babySteps(ulong thePrime, ulong _q, ulong _r2, ulong _one, ushort *h_table, ushort *h_olist, ulong *h_BJ64)\n" \
"{\n" \
"uint    j;\n" \
"ulong   resBJ, resBJ0;\n" \
"ulong   resBase = mmmNToRes(BASE, thePrime, _q, _r2);\n" \
"resBJ = resBJ0 = mmmPowmod(resBase, N_MIN, thePrime, _q, _one);\n" \
"for (j=0; j<BABY_STEPS; j++)\n" \
"{\n" \
"hashInsert(resBJ, j, h_table, h_olist, h_BJ64);\n" \
"resBJ = mmmMulmod(resBJ, resBase, thePrime, _q);\n" \
"if (resBJ == resBJ0)\n" \
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
"#if defined(USE_OPENCL)\n" \
"void collectFactor(ulong   p,\n" \
"ulong   n,\n" \
"int     c,\n" \
"volatile __global uint   *factorCount,\n" \
"__global long4  *factors)\n" \
"#else\n" \
"void collectFactor(ulong       p,\n" \
"ulong       n,\n" \
"int         c,\n" \
"volatile device atomic_int *factorCount,\n" \
"device long4      *factors)\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int old = atomic_inc(factorCount);\n" \
"#else\n" \
"int old = atomic_fetch_add_explicit(factorCount, 1, memory_order_relaxed);\n" \
"#endif\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = n;\n" \
"factors[old].y = c;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
