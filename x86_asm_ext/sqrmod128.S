/* sqrmod128.S -- (C) Geoffrey Reynolds, March 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void sqrmod128(uint64_t *R, const uint64_t *N, uint64_t inv);

   Assign R <-- R^2 (mod N) in Montgomery form, where 0 <= R < N < 2^128.
   Assumes that R is in Montgomery form, N is odd, and N[0]*inv = -1.


   void sqrmod128_proth0(uint64_t *R, const uint64_t *N, uint64_t inv);

   As sqrmod128, but assumes that N[1] < 2^63.


   void sqrmod128_proth1(uint64_t *R, const uint64_t *N);

   As sqrmod128, but assumes that N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.text
	.globl	_sqrmod128
	.globl	sqrmod128
	.p2align 4,,15

_sqrmod128:
sqrmod128:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
#endif
	mov	(%rdi), %rax
	mov	(%rsi), %r9
	mov	8(%rsi), %r10

	push	%rbp
	push	%rbx
	mov	%rdx, %rsi

	/* %rdi = R
	   %rsi = inv
	   %r9  = N[0]
	   %r10 = N[1] */

	mul	%rax
	mov	%rax, %rbx
	mov	%rdx, %rcx
	imul	%rsi, %rax
	mov	%rax, %r8

	mul	%r9
	mov	%rax, %rbp
	mov	%rdx, %r11

	mov	%r10, %rax
	mul	%r8
	add	%rax, %r11
	adc	$0, %rdx
	mov	%rdx, %r8

	mov	(%rdi), %rax
	mulq	8(%rdi)
#ifdef _WIN64
	mov	%rax, 40(%rsp)
	mov	%rdx, 48(%rsp)
#else
	mov	%rax, -8(%rsp)
	mov	%rdx, -16(%rsp)
#endif
	add	%rbp, %rbx
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%ebp, %ebp
	add	%rcx, %r11
	adc	%rdx, %r8
	adc	$0, %ebp


	add	%r11, %rax
	imul	%rsi, %rax
	mov	%rax, %rsi

	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	%r10, %rax
	mul	%rsi
	add	%rbx, %r11
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%esi, %esi
	add	%rcx, %r8
	adc	%rdx, %rbp
	adc	$0, %esi

	mov	8(%rdi), %rax
#ifdef _WIN64
	mov	40(%rsp), %rbx
	mov	48(%rsp), %rcx
#else
	mov	-8(%rsp), %rbx
	mov	-16(%rsp), %rcx
#endif

	mul	%rax
	add	%rbx, %r11
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%r8, %rcx
	adc	%rbp, %rdx
	adc	$0, %esi

0:	mov	%rcx, (%rdi)
	mov	%rdx, 8(%rdi)
	sub	%r9, %rcx
	sbb	%r10, %rdx
	sbb	$0, %esi
	jnc	0b

	pop	%rbx
	pop	%rbp
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret



	.text
	.globl	_sqrmod128_proth0
	.globl	sqrmod128_proth0
	.p2align 4,,15

_sqrmod128_proth0:
sqrmod128_proth0:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
#endif
	mov	(%rdi), %rax
	mov	(%rsi), %r9
	mov	8(%rsi), %r10
	mov	%rdx, %r8

	push	%rbp
	push	%rbx

	/* %rdi = R
	   %r8 = inv
	   %r9  = N[0]
	   %r10 = N[1] < 2^63 */

	mul	%rax
	mov	%rax, %rbx
	mov	%rdx, %rcx
	imul	%r8, %rax
	mov	%rax, %rsi

	mul	%r9
	mov	%rax, %rbp
	mov	%rdx, %r11

	mov	%rsi, %rax
	xor	%esi, %esi
	mul	%r10		/* < 2^63 */
	add	%rax, %r11
	adc	%rdx, %rsi	/* < 2^63 */

	mov	(%rdi), %rax
	mulq	8(%rdi)		/* < 2^63-1 */
	add	%rbp, %rbx
	adc	%rax, %rcx
	adc	%rdx, %rsi	/* < 2^64-1 */
	add	%rcx, %r11
	adc	$0, %rsi	/* CF=0 */
	mov	%rax, %rbx
	mov	%rdx, %rcx	/* < 2^63-1 */

	add	%r11, %rax
	imul	%r8, %rax
	mov	%rax, %r8

	mul	%r9
	xor	%ebp, %ebp
	add	%rax, %r11
	adc	%rdx, %rsi
	adc	$0, %rbp	/* < 2 */

	mov	%r10, %rax	/* < 2^63 */
	mul	%r8
	add	%rax, %rsi
	adc	%rdx, %rbp	/* < 2^63 */

	mov	8(%rdi), %rax	/* < 2^63 */
	mul	%rax
	add	%rbx, %r11
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^62 */
	add	%rsi, %rcx
	adc	%rbp, %rdx	/* CF=0 */

0:	mov	%rcx, (%rdi)
	mov	%rdx, 8(%rdi)
	sub	%r9, %rcx
	sbb	%r10, %rdx
	jnc	0b

	pop	%rbx
	pop	%rbp
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret



	.text
	.globl	_sqrmod128_proth1
	.globl	sqrmod128_proth1
	.p2align 4,,15

_sqrmod128_proth1:
sqrmod128_proth1:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
#endif
	push	%rbp
	mov	(%rdi), %rax
	mov	8(%rsi), %r10

	/* %rdi = R
	   %r10 = N[1]  */

	xor	%r9, %r9
	mul	%rax
	neg	%rax
	adc	%rdx, %r9

	mul	%r10
	xor	%r11, %r11
	add	%rax, %r9
	adc	%rdx, %r11	/* CF=0 */

	mov	(%rdi), %rax
	mulq	8(%rdi)
	xor	%ebp, %ebp
	add	%rax, %r9
	adc	%rdx, %r11
	adc	$0, %ebp	/* < 2 */
	mov	%rax, %r8
	mov	%rdx, %rsi

	add	%r9, %rax
	neg	%rax
	mov	%rax, %rcx

	mul	%r10
	add	%rcx, %r9
	adc	%rax, %r11
	adc	%rdx, %rbp	/* CF=0 */

	mov	8(%rdi), %rax
	mul	%rax
	add	%r8, %r9
	adc	%rax, %rsi
	adc	$0, %rdx	/* CF=0 */
	xor	%ecx, %ecx
	add	%r11, %rsi
	adc	%rbp, %rdx
	adc	$0, %ecx	/* < 2 */

0:	mov	%rsi, (%rdi)
	mov	%rdx, 8(%rdi)
	sub	$1, %rsi
	sbb	%r10, %rdx
	sbb	$0, %ecx
	jnc	0b

	pop	%rbp
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
