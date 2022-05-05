#ifdef __cplusplus
#ifndef _PIX_KERNEL_CL
#define _PIX_KERNEL_CL
const char *pix_kernel= \
"#if defined(USE_OPENCL)\n" \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MUL_HI          mul_hi\n" \
"void collectFactor(uint    termLength,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong2 *factors);\n" \
"#else\n" \
"#define MUL_HI          mulhi\n" \
"void collectFactor(uint        termLength,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong2     *factors);\n" \
"#endif\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"__kernel void pix_kernel(__global const  ulong  *primes,\n" \
"__global        ulong  *residuals,\n" \
"__global        uint   *terms,\n" \
"volatile __global        uint   *factorCount,\n" \
"__global        ulong2 *factors)\n" \
"{\n" \
"int   gid = get_global_id(0);\n" \
"ulong thePrime = primes[gid];\n" \
"ulong pOne = mmmOne(thePrime);\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, pOne);\n" \
"ulong multiplier = terms[0];\n" \
"ulong lengthOfEachTerm = terms[1];\n" \
"uint  termLength = terms[2];\n" \
"ulong res = (termLength == 0 ? 0 : residuals[gid]);\n" \
"ulong resMultiplier = mmmNtoRes(multiplier, thePrime, _q, _r2);\n" \
"uint dlIdx = 3;\n" \
"if (multiplier == 10)\n" \
"{\n" \
"ulong resDigit[10];\n" \
"resDigit[0] = 0;\n" \
"resDigit[1] = pOne;\n" \
"resDigit[2] = mmmAdd(resDigit[1], pOne, thePrime);\n" \
"resDigit[3] = mmmAdd(resDigit[2], pOne, thePrime);\n" \
"resDigit[4] = mmmAdd(resDigit[3], pOne, thePrime);\n" \
"resDigit[5] = mmmAdd(resDigit[4], pOne, thePrime);\n" \
"resDigit[6] = mmmAdd(resDigit[5], pOne, thePrime);\n" \
"resDigit[7] = mmmAdd(resDigit[6], pOne, thePrime);\n" \
"resDigit[8] = mmmAdd(resDigit[7], pOne, thePrime);\n" \
"resDigit[9] = mmmAdd(resDigit[8], pOne, thePrime);\n" \
"while (terms[dlIdx] != D_EOL)\n" \
"{\n" \
"res = mmmMulmod(res, resMultiplier, thePrime, _q);\n" \
"res = mmmAdd(res, resDigit[terms[dlIdx]], thePrime);\n" \
"termLength += lengthOfEachTerm;\n" \
"if (res == 0)\n" \
"collectFactor(termLength, thePrime, factorCount, factors);\n" \
"dlIdx++;\n" \
"}\n" \
"}\n" \
"else\n" \
"{\n" \
"ulong resAdder;\n" \
"while (terms[dlIdx] != D_EOL)\n" \
"{\n" \
"res = mmmMulmod(res, resMultiplier, thePrime, _q);\n" \
"resAdder = mmmNtoRes(terms[dlIdx], thePrime, _q, _r2);\n" \
"res = mmmAdd(res, resAdder, thePrime);\n" \
"termLength += lengthOfEachTerm;\n" \
"if (res == 0)\n" \
"collectFactor(termLength, thePrime, factorCount, factors);\n" \
"dlIdx++;\n" \
"}\n" \
"}\n" \
"residuals[gid] = res;\n" \
"}\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
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
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q)\n" \
"{\n" \
"ulong lo, hi;\n" \
"#ifdef __NV_CL_C_VERSION\n" \
"const uint a0 = (uint)(a), a1 = (uint)(a >> 32);\n" \
"const uint b0 = (uint)(b), b1 = (uint)(b >> 32);\n" \
"uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;\n" \
"asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a0), \"r\" (b1), \"r\" (c1));\n" \
"asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c2) : \"r\" (a0), \"r\" (b1));\n" \
"asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b1), \"r\" (c2));\n" \
"asm volatile (\"madc.hi.u32 %0, %1, %2, 0;\" : \"=r\" (c3) : \"r\" (a1), \"r\" (b1));\n" \
"asm volatile (\"mad.lo.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c1) : \"r\" (a1), \"r\" (b0), \"r\" (c1));\n" \
"asm volatile (\"madc.hi.cc.u32 %0, %1, %2, %3;\" : \"=r\" (c2) : \"r\" (a1), \"r\" (b0), \"r\" (c2));\n" \
"asm volatile (\"addc.u32 %0, %1, 0;\" : \"=r\" (c3) : \"r\" (c3));\n" \
"lo = upsample(c1, c0); hi = upsample(c3, c2);\n" \
"#else\n" \
"lo = a * b; hi = MUL_HI(a, b);\n" \
"#endif\n" \
"ulong m = lo * _q;\n" \
"ulong mp = MUL_HI(m, _p);\n" \
"long r = (long)(hi - mp);\n" \
"return (r < 0) ? r + _p : r;\n" \
"}\n" \
"#if defined(USE_OPENCL)\n" \
"void collectFactor(uint    termLength,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong2 *factors)\n" \
"#else\n" \
"void collectFactor(uint        termLength,\n" \
"ulong       p,\n" \
"volatile device atomic_int *factorCount,\n" \
"device ulong2     *factors)\n" \
"#endif\n" \
"{\n" \
"#if defined(USE_OPENCL)\n" \
"int old = atomic_inc(factorCount);\n" \
"#else\n" \
"int old = atomic_fetch_add_explicit(factorCount, 1, memory_order_relaxed);\n" \
"#endif\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = termLength;\n" \
"factors[old].y = p;\n" \
"}\n" \
;
#endif
#endif
