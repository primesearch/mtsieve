#ifdef __cplusplus
#ifndef _PRIMORIAL_KERNEL_CL
#define _PRIMORIAL_KERNEL_CL
const char *primorial_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"void collectFactor(uint    primorial,\n" \
"int     c,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorsCount,\n" \
"__global long4  *factors);\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmSub(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q);\n" \
"__kernel void primorial_kernel(__global const  ulong  *primes,\n" \
"__global const  ushort *primorialGaps,\n" \
"__global const  ulong  *params,\n" \
"__global        ulong2 *rems,\n" \
"volatile __global        uint   *factorsCount,\n" \
"__global        long4  *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"int    steps;\n" \
"int    idx = 0;\n" \
"int    maxIdxFirstLoop = D_MIN_IDX;\n" \
"ulong  resGaps[D_BIGGEST_GAP+1];\n" \
"ulong  currentIdx = params[0];\n" \
"ulong  thePrime = primes[gid];\n" \
"if (thePrime == 0) return;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong pOne = mmmOne(thePrime);\n" \
"ulong mOne = mmmSub(0, pOne, thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, pOne);\n" \
"ulong ri, rf;\n" \
"resGaps[0] = 0;\n" \
"resGaps[2] = mmmNtoRes(2, thePrime, _q, _r2);\n" \
"for (idx=4; idx<=D_BIGGEST_GAP; idx+=2)\n" \
"resGaps[idx] = mmmAdd(resGaps[idx-2], resGaps[2], thePrime);\n" \
"if (currentIdx == 0)\n" \
"{\n" \
"ri = mmmNtoRes(FIRST_PRIMORIAL_PRIME, thePrime, _q, _r2);\n" \
"rf = mmmNtoRes(FIRST_PRIMORIAL, thePrime, _q, _r2);\n" \
"}\n" \
"else\n" \
"{\n" \
"ri = rems[gid].x;\n" \
"rf = rems[gid].y;\n" \
"}\n" \
"steps = 0;\n" \
"idx = currentIdx;\n" \
"for (; idx<maxIdxFirstLoop; idx++, steps++)\n" \
"{\n" \
"if (steps >= D_MAX_STEPS)\n" \
"break;\n" \
"ri = mmmAdd(ri, resGaps[primorialGaps[idx]], thePrime);\n" \
"rf = mmmMulmod(rf, ri, thePrime, _q);\n" \
"}\n" \
"while (primorialGaps[idx] > 0)\n" \
"{\n" \
"if (steps >= D_MAX_STEPS)\n" \
"break;\n" \
"ri = mmmAdd(ri, resGaps[primorialGaps[idx]], thePrime);\n" \
"rf = mmmMulmod(rf, ri, thePrime, _q);\n" \
"if (rf == pOne)\n" \
"collectFactor(mmmResToN(ri, thePrime, _q), -1, thePrime, factorsCount, factors);\n" \
"if (rf == mOne)\n" \
"collectFactor(mmmResToN(ri, thePrime, _q), +1, thePrime, factorsCount, factors);\n" \
"idx++;\n" \
"steps++;\n" \
"}\n" \
"rems[gid].x = ri;\n" \
"rems[gid].y = rf;\n" \
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
"lo = a * b; hi = mul_hi(a, b);\n" \
"#endif\n" \
"ulong m = lo * _q;\n" \
"ulong mp = mul_hi(m, _p);\n" \
"long r = (long)(hi - mp);\n" \
"return (r < 0) ? r + _p : r;\n" \
"}\n" \
"ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
"}\n" \
"ulong mmmResToN(ulong res, ulong _p, ulong _q)\n" \
"{\n" \
"return mmmMulmod(res, 1, _p, _q);\n" \
"}\n" \
"void collectFactor(uint    primorial,\n" \
"int     c,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorsCount,\n" \
"__global long4  *factors)\n" \
"{\n" \
"if (primorial < D_MIN_PRIMORIAL)\n" \
"return;\n" \
"int old = atomic_inc(factorsCount);\n" \
"if (old >= D_MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = primorial;\n" \
"factors[old].y = c;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
