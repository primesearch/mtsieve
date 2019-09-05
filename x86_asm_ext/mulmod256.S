/* mulmod256.S -- (C) Geoffrey Reynolds, June 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void mulmod256(uint64_t *R, const uint64_t *A, const uint64_t *B,
                  const uint64_t *N, uint64_t inv);

   Assign R <-- A*B (mod N) in Montgomery form, where 0 <= A,B < N < 2^256.
   Assumes that A,B are in Montgomery form, N is odd, and N[0]*inv = -1.


   void mulmod256_proth0(uint64_t *R, const uint64_t *A, const uint64_t *B,
                         const uint64_t *N);

   As above, but assumes that N[3] < 2^63, N[1] = 0, N[0] = 1.


   void mulmod256_proth1(uint64_t *R, const uint64_t *A, const uint64_t *B,
                         const uint64_t *N);

   As above, but assumes that N[2] = N[1] = 0, N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#if defined(_WIN64) || defined(__APPLE__)
# define FUN_NAME _mulmod256
# define FUN_NAME0 _mulmod256_proth0
# define FUN_NAME1 _mulmod256_proth1
#else
# define FUN_NAME mulmod256
# define FUN_NAME0 mulmod256_proth0
# define FUN_NAME1 mulmod256_proth1
#endif

	.globl	_mulmod256
	.globl	mulmod256
	.p2align 4,,15

_mulmod256:
mulmod256:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
	mov	64(%rsp), %r8
#endif
	mov	%rdx, %r11
	mov	%rcx, %r10

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15
	push	%r8

	/* (%rsp) = inv
	   %r10 = N
	   %r11 = B
	   %rsi = A
	   %rdi = R */

	mov	(%rsi), %r9
	imul	(%r11), %r9
	imul	%r8, %r9

	mov	(%r10), %rax
	mul	%r9
	mov	%rax, %r14
	mov	%rdx, %r15

	mov	8(%r10), %rax
	mul	%r9
	xor	%r8, %r8
	add	%rax, %r15
	adc	%rdx, %r8

	mov	16(%r10), %rax
	mul	%r9
	xor	%rbp, %rbp
	add	%rax, %r8
	adc	%rdx, %rbp

	mov	24(%r10), %rax
	mul	%r9
	xor	%r12, %r12
	add	%rax, %rbp
	adc	%rdx, %r12

	mov	(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r14
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r15
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rcx, %rbp
	adc	%rdx, %r12
	adc	$0, %r13


	mov	8(%rsi), %r9
	imul	(%r11), %r9
	add	%r15, %r9
	imul	(%rsp), %r9

	mov	(%r10), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r10), %rax
	mul	%r9
	add	%rbx, %r15
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r10), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r10), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx
	add	%rcx, %r12
	adc	$0, %rbx

	xor	%r14, %r14
	add	%rbx, %r13
	adc	$0, %r14
	mov	8(%rsi), %r9

	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r15
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r12
	adc	%rdx, %r13
	adc	$0, %r14


	mov	16(%rsi), %r9
	imul	(%r11), %r9
	add	%r8, %r9
	imul	(%rsp), %r9

	mov	(%r10), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r10), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r10), %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r10), %rax
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx
	add	%rcx, %r13
	adc	$0, %rbx

	xor	%r15, %r15
	add	%rbx, %r14
	adc	$0, %r15
	mov	16(%rsi), %r9

	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r13
	adc	%rdx, %r14
	adc	$0, %r15


	mov	24(%rsi), %r9
	imul	(%r11), %r9
	add	%rbp, %r9
	imul	(%rsp), %r9

	mov	(%r10), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r10), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r10), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r10), %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx
	add	%rcx, %r14
	adc	$0, %rbx

	xor	%r8, %r8	/* carry */
	add	%rbx, %r15
	adc	$0, %r8
	mov	24(%rsi), %r9

	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r14
	adc	%rdx, %r15
	adc	$0, %r8

0:
	mov	%r12, (%rdi)
	mov	%r13, 8(%rdi)
	mov	%r14, 16(%rdi)
	mov	%r15, 24(%rdi)
	sub	(%r10), %r12
	sbb	8(%r10), %r13
	sbb	16(%r10), %r14
	sbb	24(%r10), %r15
	sbb	$0, %r8
	jnc	0b

	pop	%r8
	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
   


	.text
	.globl	_mulmod256_proth0
	.globl	mulmod256_proth0
	.p2align 4,,15

_mulmod256_proth0:
mulmod256_proth0:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
#endif
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	mov	%rdx, %r11
	mov	%rcx, %r10

	/* %r10 = N
	   %r11 = B
	   %rsi = A
	   %rdi = R */

	mov	(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rdx, %r15

	neg	%rax
	mov	%rax, %r14

	mov	16(%r10), %rax
	mul	%r14
	mov	%rax, %r8
	mov	%rdx, %rbp

	mov	24(%r10), %rax	/* < 2^63 */
	mul	%r14
	xor	%r12, %r12
	add	%rax, %rbp
	adc	%rdx, %r12	/* < 2^63 */

	mov	8(%r11), %rax
	mul	%r9
	neg	%r14
	adc	%rax, %r15
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	xor	%rcx, %rcx
	add	%rax, %rbx
	adc	%rdx, %rcx

	mov	24(%r11), %rax	/* < 2^63-1 */
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rcx, %rbp
	adc	%rdx, %r12	/* CF=0 */


	mov	8(%rsi), %r9
	imul	(%r11), %r9
	add	%r15, %r9
	neg	%r9

	mov	16(%r10), %rax
	mul	%r9
	xor	%rcx, %rcx
	add	%r9, %r15
	adc	$0, %r8
	adc	%rax, %rbp
	adc	%rdx, %rcx

	mov	24(%r10), %rax	/* < 2^63 */
	mul	%r9
	xor	%rbx, %rbx
	add	%rax, %rcx
	adc	%rdx, %rbx	/* < 2^63 */
	xor	%r13, %r13
	add	%rcx, %r12
	adc	%rbx, %r13	/* <= 2^63, so CF=0 */

	mov	8(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r15
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax	/* < 2^63-1 */
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rcx, %r12
	adc	%rdx, %r13	/* CF=0 */


	mov	16(%rsi), %r9
	imul	(%r11), %r9
	add	%r8, %r9
	neg	%r9

	mov	16(%r10), %rax
	mul	%r9
	xor	%rbx, %rbx
	xor	%rcx, %rcx
	add	%r9, %r8
	adc	$0, %rbp
	adc	%rax, %rbx
	adc	%rdx, %rcx

	mov	24(%r10), %rax	/* < 2^63 */
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %r14	/* < 2^63-1 */
	adc	$0, %r14	/* < 2^63 */
	add	%rcx, %r13
	adc	$0, %r14	/* <= 2^63 */

	mov	16(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax	/* < 2^63-1 */
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rcx, %r13
	adc	%rdx, %r14	/* CF=0 */


	mov	24(%rsi), %r9
	imul	(%r11), %r9
	add	%rbp, %r9
	neg	%r9

	mov	16(%r10), %rax
	mul	%r9
	xor	%rbx, %rbx
	xor	%rcx, %rcx
	add	%r9, %rbp
	adc	$0, %r12
	adc	%rax, %rbx
	adc	%rdx, %rcx

	mov	24(%r10), %rax	/* < 2^63 */
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	mov	%rdx, %r15	/* < 2^63-1 */
	adc	$0, %r15	/* < 2^63 */
	add	%rcx, %r14
	adc	$0, %r15	/* <= 2^63 */

	mov	24(%rsi), %r9	/* < 2^63-1 */
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx	/* < 2^63-2 */

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx	/* < 2^63-2 */
	adc	$0, %rbx	/* < 2^63-1 */

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx	/* < 2^63-2 */
	adc	$0, %rcx	/* < 2^63-1 */

	mov	24(%r11), %rax	/* < 2^63-1 */
	mul	%r9		/* < 2^62-1 */
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^62 */
	add	%rcx, %r14
	adc	%rdx, %r15	/* CF=0 */

0:
	mov	%r12, (%rdi)
	mov	%r13, 8(%rdi)
	mov	%r14, 16(%rdi)
	mov	%r15, 24(%rdi)
	sub	$1, %r12
	sbb	$0, %r13
	sbb	16(%r10), %r14
	sbb	24(%r10), %r15
	jnc	0b

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret


	.text
	.globl	_mulmod256_proth1
	.globl	mulmod256_proth1
	.p2align 4,,15

_mulmod256_proth1:
mulmod256_proth1:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
#endif
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	mov	%rdx, %r11
	mov	24(%rcx), %r10

	/* %r10 = N[3]
	   %r11 = B
	   %rsi = A
	   %rdi = R */

	mov	(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %r15

	neg	%rax

	mul	%r10
	mov	%rax, %rbp
	mov	%rdx, %r12

	mov	8(%r11), %rax
	mul	%r9
	xor	%r8, %r8
	neg	%rbx
	adc	%rax, %r15
	adc	%rdx, %r8

	mov	16(%r11), %rax
	mul	%r9
	xor	%rcx, %rcx
	add	%rax, %r8
	adc	%rdx, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rax, %rcx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rcx, %rbp
	adc	%rdx, %r12
	adc	$0, %r13


	mov	8(%rsi), %r9
	imul	(%r11), %r9
	add	%r15, %r9
	neg	%r9

	mov	%r10, %rax
	mul	%r9
	xor	%r14, %r14
	add	%r9, %r15
	adc	$0, %r8
	adc	$0, %rbp
	adc	%rax, %r12
	adc	%rdx, %r13
	adc	$0, %r14

	mov	8(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r15
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r12
	adc	%rdx, %r13
	adc	$0, %r14


	mov	16(%rsi), %r9
	imul	(%r11), %r9
	add	%r8, %r9
	neg	%r9

	mov	%r10, %rax
	mul	%r9
	xor	%r15, %r15
	add	%r9, %r8
	adc	$0, %rbp
	adc	$0, %r12
	adc	%rax, %r13
	adc	%rdx, %r14
	adc	$0, %r15

	mov	16(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r13
	adc	%rdx, %r14
	adc	$0, %r15


	mov	24(%rsi), %r9
	imul	(%r11), %r9
	add	%rbp, %r9
	neg	%r9

	mov	%r10, %rax
	mul	%r9

	xor	%r8, %r8	/* carry */
	add	%r9, %rbp
	adc	$0, %r12
	adc	$0, %r13
	adc	%rax, %r14
	adc	%rdx, %r15
	adc	$0, %r8

	mov	24(%rsi), %r9
	mov	(%r11), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r11), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r11), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%r11), %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r14
	adc	%rdx, %r15
	adc	$0, %r8

0:
	mov	%r12, (%rdi)
	mov	%r13, 8(%rdi)
	mov	%r14, 16(%rdi)
	mov	%r15, 24(%rdi)
	sub	$1, %r12
	sbb	$0, %r13
	sbb	$0, %r14
	sbb	%r10, %r15
	sbb	$0, %r8
	jnc	0b

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
