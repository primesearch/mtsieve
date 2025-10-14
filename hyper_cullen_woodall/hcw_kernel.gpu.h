#ifdef __cplusplus
#ifndef _HCW_KERNEL_CL
#define _HCW_KERNEL_CL
const char *hcw_kernel= \
"#define MAX_POWERS 20\n" \
"#if defined(USE_OPENCL)\n" \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MUL_HI          mul_hi\n" \
"void collectFactor(uint    b,\n" \
"uint    n,\n" \
"int     sign,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"#else\n" \
"void collectFactor(uint        b,\n" \
"uint        n,\n" \
"int         sign,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong4     *factors);\n" \
"#endif\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);\n" \
"__kernel void hcw_kernel(__global const ulong  *primes,\n" \
"__global const uint   *terms,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"uint   gid = get_global_id(0);\n" \
"uint   termIndex;\n" \
"uint   currentB;\n" \
"uint   currentN, previousN;\n" \
"uint   pIndex;\n" \
"ulong  thePrime = primes[gid];\n" \
"ulong  magicNumber;\n" \
"ulong  magicShift;\n" \
"ulong  resB, resN;\n" \
"ulong  resBexpN, resTemp;\n" \
"ulong  resNexpB;\n" \
"ulong  resPowers[MAX_POWERS];\n" \
"if (thePrime == 0)\n" \
"return;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong _zero = mmmNToRes(0, thePrime, _q, _r2);\n" \
"#ifdef IS_PLUS\n" \
"ulong _pm1 = mmmSub(_zero, _one, thePrime);\n" \
"#endif\n" \
"termIndex = 0;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"currentB = terms[termIndex];\n" \
"resB = mmmNToRes(currentB, thePrime, _q, _r2);\n" \
"termIndex++;\n" \
"if (terms[termIndex] > 0)\n" \
"{\n" \
"currentN = terms[termIndex];\n" \
"resN = mmmNToRes(currentN, thePrime, _q, _r2);\n" \
"resPowers[0] = _one;\n" \
"resPowers[1] = resB;\n" \
"for (pIndex=2; pIndex<=MAX_POWERS; pIndex++)\n" \
"resPowers[pIndex] = mmmMulmod(resPowers[pIndex-1], resB, thePrime, _q);\n" \
"resBexpN = mmmPowmod(resB, currentN, thePrime, _q, _one);\n" \
"resNexpB = mmmPowmod(resN, currentB, thePrime, _q, _one);\n" \
"ulong res = mmmMulmod(resBexpN, resNexpB, thePrime, _q);\n" \
"#ifdef IS_MINUS\n" \
"if (res == _one)\n" \
"collectFactor(currentB, currentN, -1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"#ifdef IS_PLUS\n" \
"if (res == _pm1)\n" \
"collectFactor(currentB, currentN, +1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"previousN = currentN;\n" \
"termIndex++;\n" \
"while (terms[termIndex] > 0)\n" \
"{\n" \
"currentN = terms[termIndex];\n" \
"pIndex = currentN - previousN;\n" \
"if (pIndex > MAX_POWERS)\n" \
"resTemp = mmmPowmod(resB, currentN - previousN, thePrime, _q, _one);\n" \
"else\n" \
"resTemp = resPowers[pIndex];\n" \
"resBexpN = mmmMulmod(resBexpN, resTemp, thePrime, _q);\n" \
"resN = mmmNToRes(currentN, thePrime, _q, _r2);\n" \
"resNexpB = mmmPowmod(resN, currentB, thePrime, _q, _one);\n" \
"ulong res = mmmMulmod(resBexpN, resNexpB, thePrime, _q);\n" \
"#ifdef IS_MINUS\n" \
"if (res == _one)\n" \
"collectFactor(currentB, currentN, -1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"#ifdef IS_PLUS\n" \
"if (res == _pm1)\n" \
"collectFactor(currentB, currentN, +1, thePrime, factorCount, factors);\n" \
"#endif\n" \
"termIndex++;\n" \
"previousN = currentN;\n" \
"}\n" \
"}\n" \
"termIndex++;\n" \
"};\n" \
"}\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q)\n" \
"{\n" \
"return mmmMulmod(res, 1, _p, _q);\n" \
"}\n" \
"ulong mmmInvert(ulong _p)\n" \
"{\n" \
"ulong p_inv = 1;\n" \
"ulong prev = 0;\n" \
"while (p_inv != prev)\n" \
"{\n" \
"prev = p_inv;\n" \
"p_inv *= (2 - _p * p_inv);\n" \
"}\n" \
"return p_inv;\n" \
"}\n" \
"ulong mmmOne(ulong _p)\n" \
"{\n" \
"return ((-_p) % _p);\n" \
"}\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one)\n" \
"{\n" \
"ulong t = mmmAdd(_one, _one, _p);\n" \
"t = mmmAdd(t, t, _p);\n" \
"for (size_t i=0; i<5; i++)\n" \
"t = mmmMulmod(t, t, _p, _q);\n" \
"return t;\n" \
"}\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a >= _p - b) ? _p : 0;\n" \
"return a + b - c;\n" \
"}\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p)\n" \
"{\n" \
"ulong c = (a < b) ? _p : 0;\n" \
"return a - b + c;\n" \
"}\n" \
"ulong mmmNToRes(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
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
"#if defined(USE_OPENCL)\n" \
"void collectFactor(uint    b,\n" \
"uint    n,\n" \
"int     sign,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"#else\n" \
"void collectFactor(uint        b,\n" \
"uint        n,\n" \
"int         sign,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device long4      *factors)\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int old = atomic_inc(factorCount);\n" \
"#else\n" \
"int old = atomic_fetch_add_explicit(factorCount, 1, memory_order_relaxed);\n" \
"#endif\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = b;\n" \
"factors[old].y = n;\n" \
"if (sign > 2)\n" \
"factors[old].z = sign;\n" \
"else\n" \
"factors[old].z = (sign > 0 ? +1 : +2);\n" \
"factors[old].w = p;\n" \
"}\n" \
;
#endif
#endif
