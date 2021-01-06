#ifdef __cplusplus
#ifndef _GFN_KERNEL_CL
#define _GFN_KERNEL_CL
const char *gfn_kernel= \
"#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable\n" \
"void collectFactor(ulong   k,\n" \
"uint    n,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"ulong mmmInvert(ulong _p);\n" \
"ulong mmmOne(ulong _p);\n" \
"ulong mmmR2(ulong _p, ulong _q, ulong _one);\n" \
"ulong mmmAdd(ulong a, ulong b, ulong _p);\n" \
"ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);\n" \
"ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2);\n" \
"ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one, ulong _r2);\n" \
"__kernel void gfn_kernel(__global const ulong  *primes,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"int    gid = get_global_id(0);\n" \
"uint   n;\n" \
"ulong  p = primes[gid];\n" \
"ulong  k = (1+p) >> 1;\n" \
"ulong _q = mmmInvert(p);\n" \
"ulong _one = mmmOne(p);\n" \
"ulong _r2 = mmmR2(p, _q, _one);\n" \
"ulong mpK = mmmN(k, p, _q, _r2);\n" \
"ulong mpRes = mmmPowmod(mpK, N_MIN, p, _q, _one, _r2);\n" \
"ulong rem = mmmMulmod(mpRes, 1, p, _q);\n" \
"k = p - rem;\n" \
"for (n=N_MIN; n<=N_MAX; n++)\n" \
"{\n" \
"if (k & 1)\n" \
"{\n" \
"if (k <= K_MAX && k >= K_MIN)\n" \
"collectFactor(k, n, p, factorCount, factors);\n" \
"k += p;\n" \
"}\n" \
"k >>= 1;\n" \
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
"ulong mmmN(ulong n, ulong _p, ulong _q, ulong _r2)\n" \
"{\n" \
"return mmmMulmod(n, _r2, _p, _q);\n" \
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
"void collectFactor(ulong   k,\n" \
"uint    n,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"{\n" \
"int old = atomic_inc(factorCount);\n" \
"if (old >= MAX_FACTORS)\n" \
"return;\n" \
"factors[old].x = k;\n" \
"factors[old].y = n;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
