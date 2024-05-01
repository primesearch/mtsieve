#ifdef __cplusplus
#ifndef _LIFCHITZ_KERNEL_CL
#define _LIFCHITZ_KERNEL_CL
const char *lifchitz_kernel= \
"#if defined(USE_OPENCL)\n" \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MUL_HI          mul_hi\n" \
"void collectFactor(ulong   idx,\n" \
"int     sign,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"#else\n" \
"void collectFactor(ulong       idx,\n" \
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
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"__kernel void lifchitz_kernel(__global const ulong  *primes,\n" \
"__global const uint   *xBases,\n" \
"__global const uint   *yBases,\n" \
"__global const uint4  *terms,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"uint   gid = get_global_id(0);\n" \
"ulong  thePrime = primes[gid];\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong resBase, res = 0;\n" \
"ulong xRemainders[D_MAX_X_BASES];\n" \
"ulong yRemainders[D_MAX_Y_BASES];\n" \
"uint  x, minX = xBases[0];\n" \
"uint  y, minY = yBases[0];\n" \
"resBase = mmmNtoRes(minX, thePrime, _q, _r2);\n" \
"for (uint xIdx=0; xIdx<D_MAX_X_BASES; xIdx++)\n" \
"{\n" \
"if (xBases[xIdx] > 0)\n" \
"{\n" \
"x = xBases[xIdx];\n" \
"res = mmmPowmod(resBase, x, thePrime, _q, _one, _r2);\n" \
"xRemainders[x-minX] = mmmResToN(res, thePrime, _q);\n" \
"}\n" \
"resBase = mmmAdd(resBase, _one, thePrime);\n" \
"}\n" \
"resBase = mmmNtoRes(minY, thePrime, _q, _r2);\n" \
"for (uint yIdx=0; yIdx<D_MAX_Y_BASES; yIdx++)\n" \
"{\n" \
"if (yBases[yIdx] > 0)\n" \
"{\n" \
"y = yBases[yIdx];\n" \
"res = mmmPowmod(resBase, y, thePrime, _q, _one, _r2);\n" \
"yRemainders[y-minY] = mmmResToN(res, thePrime, _q);\n" \
"}\n" \
"resBase = mmmAdd(resBase, _one, thePrime);\n" \
"}\n" \
"for (uint idx=0; terms[idx].x>0; idx++)\n" \
"{\n" \
"ulong xPowX = xRemainders[terms[idx].x - minX];\n" \
"ulong yPowY = yRemainders[terms[idx].y - minY];\n" \
"if (terms[idx].z == 1)\n" \
"{\n" \
"if (xPowX + yPowY == thePrime)\n" \
"collectFactor(idx, +1, thePrime, factorCount, factors);\n" \
"}\n" \
"else\n" \
"{\n" \
"if (xPowX == yPowY)\n" \
"collectFactor(idx, -1, thePrime, factorCount, factors);\n" \
"}\n" \
"}\n" \
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
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
"}\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q)\n" \
"{\n" \
"return mmmMulmod(res, 1, _p, _q);\n" \
"}\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2)\n" \
"{\n" \
"ulong x = resB;\n" \
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
"#if defined(USE_OPENCL)\n" \
"void collectFactor(ulong   idx,\n" \
"int     sign,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"#else\n" \
"void collectFactor(ulong       idx,\n" \
"int         sign,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong      *factors)\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int old = atomic_inc(factorCount);\n" \
"#else\n" \
"int old = atomic_fetch_add_explicit(factorCount, 1, memory_order_relaxed);\n" \
"#endif\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = idx;\n" \
"factors[old].y = (sign == 1 ? 1 : 2);\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
