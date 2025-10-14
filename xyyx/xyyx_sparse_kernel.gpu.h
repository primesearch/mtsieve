#ifdef __cplusplus
#ifndef _XYYX_SPARSE_KERNEL_CL
#define _XYYX_SPARSE_KERNEL_CL
const char *xyyx_sparse_kernel= \
"#define MAX_POWERS 50\n" \
"#if defined(USE_OPENCL)\n" \
"#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable\n" \
"#define MUL_HI          mul_hi\n" \
"void collectFactor(uint    x,\n" \
"uint    y,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors);\n" \
"#else\n" \
"#define MUL_HI          mulhi\n" \
"void collectFactor(uint        x,\n" \
"uint        y,\n" \
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
"ulong mmmPowmod(ulong resbase, ulong exp, ulong _p, ulong _q, ulong _one);\n" \
"__kernel void xyyx_sparse_kernel(__global const ulong  *primes,\n" \
"__global const uint2  *terms,\n" \
"volatile __global       uint   *factorCount,\n" \
"__global       ulong4 *factors)\n" \
"{\n" \
"uint   gid = get_global_id(0);\n" \
"uint   termIndex;\n" \
"ulong  thePrime = primes[gid];\n" \
"if (thePrime == 0)\n" \
"return;\n" \
"ulong _q = mmmInvert(thePrime);\n" \
"ulong _one = mmmOne(thePrime);\n" \
"ulong _r2 = mmmR2(thePrime, _q, _one);\n" \
"ulong _zero = mmmNToRes(0, thePrime, _q, _r2);\n" \
"ulong resBaseX, resBaseY;\n" \
"ulong resXexpY, resYexpX;\n" \
"for (uint idx=0; ; idx++)\n" \
"{\n" \
"if (terms[idx].x == 0)\n" \
"break;\n" \
"resBaseX = mmmNToRes(terms[idx].x, thePrime, _q, _r2);\n" \
"resBaseY = mmmNToRes(terms[idx].y, thePrime, _q, _r2);\n" \
"resXexpY = mmmPowmod(resBaseX, terms[idx].y, thePrime, _q, _one);\n" \
"resYexpX = mmmPowmod(resBaseY, terms[idx].x, thePrime, _q, _one);\n" \
"#ifdef IS_MINUS\n" \
"if (resXexpY == resYexpX)\n" \
"collectFactor(terms[idx].x, terms[idx].y, thePrime, factorCount, factors);\n" \
"#else\n" \
"resYexpX = mmmAdd(resXexpY, resYexpX, thePrime);\n" \
"if (resYexpX == _zero)\n" \
"collectFactor(terms[idx].x, terms[idx].y, thePrime, factorCount, factors);\n" \
"#endif\n" \
"};\n" \
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
"void collectFactor(uint    x,\n" \
"uint    y,\n" \
"ulong   p,\n" \
"volatile __global uint   *factorCount,\n" \
"__global ulong4 *factors)\n" \
"#else\n" \
"void collectFactor(uint        x,\n" \
"uint        y,\n" \
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
"factors[old].x = x;\n" \
"factors[old].y = y;\n" \
"factors[old].z = p;\n" \
"}\n" \
;
#endif
#endif
