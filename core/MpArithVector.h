/* MpArith.h -- (C) Mark Rodenkirch, December 2020

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
*/

#define		VECTOR_SIZE		4		// must be a power of two

// Arithmetic on vectors: hide the latency of the MUL instruction.

// Montgomery form: if 0 <= a < p then r is 2^64 * a mod p
template <size_t N>
class MpRes
{
private:
	uint64_t _r[N];

public:
	uint64_t operator [](const size_t i) const { return _r[i]; }
	uint64_t & operator [](const size_t i) { return _r[i]; }
};

// Montgomery modular arithmetic in Z/pZ
template <size_t N>
class MpArith
{
private:
	uint64_t _p[N], _q[N];
	MpRes<N> _one;		// 2^64 mod p
	MpRes<N> _r2;		// (2^64)^2 mod p

private:
	// p * p_inv = 1 (mod 2^64) (Newton's method)
	constexpr uint64_t invert(const uint64_t p)
	{
		uint64_t p_inv = 1, prev = 0;
		while (p_inv != prev) { prev = p_inv; p_inv *= 2 - p * p_inv; }
		return p_inv;
	}

	// The Montgomery REDC algorithm
	constexpr uint64_t REDC(const __uint128_t t, const uint64_t p, const uint64_t q) const
	{
		const uint64_t m = uint64_t(t) * q;
		const int64_t r = int64_t((t >> 64) - uint64_t((m * __uint128_t(p)) >> 64));
		return (r < 0) ? uint64_t(r + p) : uint64_t(r);
	}

public:
	MpArith(const uint64_t * const p)
	{
		for (size_t k = 0; k < N; ++k)
		{
			const uint64_t p_k = p[k];
			_p[k] = p_k;
			_q[k] = invert(p_k);
			_one[k] = (-p_k) % p_k;
		}
		
		MpRes<N> t = add(_one, _one); t = add(t, t);	// 4
		for (size_t i = 0; i < 5; ++i) t = mul(t, t);	// 4^{2^5} = 2^64
		_r2 = t;
	}

	// Convert n to Montgomery representation
	MpRes<N> toMp(const uint64_t n) const
	{
		// n * (2^64)^2 = (n * 2^64) * (1 * 2^64)
		MpRes<N> r;
		for (size_t k = 0; k < N; ++k) r[k] = n;
		return mul(r, _r2);
	}

	static MpRes<N> zero()
	{
		MpRes<N> r;
		for (size_t k = 0; k < N; ++k) r[k] = 0;
		return r;
	}

	MpRes<N> one() const { return _one; }	// Montgomery form of 1

	static bool at_least_one_is_equal(const MpRes<N> & a, const MpRes<N> & b)
	{
		bool is_equal = false;
		for (size_t k = 0; k < N; ++k) is_equal |= (a[k] == b[k]);
		return is_equal;
	}

	MpRes<N> add(const MpRes<N> & a, const MpRes<N> & b) const
	{
		MpRes<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			const uint64_t c = (a[k] >= _p[k] - b[k]) ? _p[k] : 0;
			r[k] = a[k] + b[k] - c;
		}
		return r;
	}

	MpRes<N> sub(const MpRes<N> & a, const MpRes<N> & b) const
	{
		MpRes<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			const uint64_t c = (a[k] < b[k]) ? _p[k] : 0;
			r[k] = a[k] - b[k] + c;
		}
		return r;
	}

	MpRes<N> mul(const MpRes<N> & a, const MpRes<N> & b) const
	{
		MpRes<N> r;
		for (size_t k = 0; k < N; ++k)
		{
			r[k] = REDC(a[k] * __uint128_t(b[k]), _p[k], _q[k]);
		}
		return r;
	}
};


typedef MpRes<VECTOR_SIZE> MpResVec;
typedef MpArith<VECTOR_SIZE> MpArithVec;
