/* sqrmod256.S -- (C) Geoffrey Reynolds, June 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void sqrmod256(uint64_t *R, const uint64_t *N, uint64_t inv);

   Assign R <-- R^2 (mod N) in Montgomery form, where 0 <= R < N < 2^256.
   Assumes that R is in Montgomery form, N is odd, and N[0]*inv = -1.


   void sqrmod256_proth0(uint64_t *R, const uint64_t *N);

   As sqrmod256, but assumes that N[3] < 2^63, N[1] = 0, N[0] = 1.


   void sqrmod256_proth1(uint64_t *R, const uint64_t *N);

   As sqrmod256, but assumes that N[2] = N[1] = 0, N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.text
	.globl	_sqrmod256
	.globl	sqrmod256
	.p2align 4,,15

_sqrmod256:
sqrmod256:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
#endif
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%rdx

	/* (%rsp) = inv
	   %rsi = N
	   %rdi = R  */

	mov	(%rdi), %r9
	imul	%r9, %r9
	imul	%rdx, %r9

	mov	(%rsi), %rax
	mul	%r9
	mov	%rax, %r10
	mov	%rdx, %r11

	mov	8(%rsi), %rax
	mul	%r9
	xor	%r8, %r8
	add	%rax, %r11
	adc	%rdx, %r8

	mov	16(%rsi), %rax
	mul	%r9
	xor	%rbp, %rbp
	add	%rax, %r8
	adc	%rdx, %rbp

	mov	24(%rsi), %rax
	mul	%r9
	xor	%r12, %r12
	add	%rax, %rbp
	adc	%rdx, %r12

	mov	(%rdi), %r9
	mov	%r9, %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %r10
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r11
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rcx, %rbp
	adc	%rdx, %r12
	adc	$0, %r13


	mov	8(%rdi), %r9
	imul	(%rdi), %r9
	add	%r11, %r9
	imul	(%rsp), %r9

	mov	(%rsi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rsi), %rax
	mul	%r9
	add	%rbx, %r11
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rsi), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rsi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx
	add	%rcx, %r12
	adc	$0, %rbx

	xor	%r10, %r10
	add	%rbx, %r13
	adc	$0, %r10
	mov	8(%rdi), %r9

	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	%r9, %rax
	mul	%r9
	add	%rbx, %r11
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r12
	adc	%rdx, %r13
	adc	$0, %r10


	mov	16(%rdi), %r9
	imul	(%rdi), %r9
	add	%r8, %r9
	imul	(%rsp), %r9

	mov	(%rsi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rsi), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rsi), %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rsi), %rax
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx
	add	%rcx, %r13
	adc	$0, %rbx

	xor	%r11, %r11
	add	%rbx, %r10
	adc	$0, %r11
	mov	16(%rdi), %r9

	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	%r9, %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r13
	adc	%rdx, %r10
	adc	$0, %r11


	mov	24(%rdi), %r9
	imul	(%rdi), %r9
	add	%rbp, %r9
	imul	(%rsp), %r9

	mov	(%rsi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rsi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rsi), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rsi), %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx
	add	%rcx, %r10
	adc	$0, %rbx

	xor	%r8, %r8	/* carry */
	add	%rbx, %r11
	adc	$0, %r8
	mov	24(%rdi), %r9

	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	%r9, %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r10
	adc	%rdx, %r11
	adc	$0, %r8

0:
	mov	%r12, (%rdi)
	mov	%r13, 8(%rdi)
	mov	%r10, 16(%rdi)
	mov	%r11, 24(%rdi)
	sub	(%rsi), %r12
	sbb	8(%rsi), %r13
	sbb	16(%rsi), %r10
	sbb	24(%rsi), %r11
	sbb	$0, %r8
	jnc	0b

	pop	%rdx
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
	.globl	_sqrmod256_proth0
	.globl	sqrmod256_proth0
	.align	16

_sqrmod256_proth0:
sqrmod256_proth0:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
#endif
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	/* %rsi = N
	   %rdi = R  */

	mov	(%rdi), %r9
	mov	%r9, %rax
	mul	%r9
	mov	%rdx, %r11

	neg	%rax
	mov	%rax, %r10

	mov	16(%rsi), %rax
	mul	%r10
	mov	%rax, %r8
	mov	%rdx, %rbp

	mov	24(%rsi), %rax	/* < 2^63 */
	mul	%r10
	xor	%r12, %r12
	add	%rax, %rbp
	adc	%rdx, %r12	/* < 2^63, so CF=0 */

	mov	8(%rdi), %rax
	mul	%r9
	neg	%r10
	adc	%rax, %r11
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	xor	%rcx, %rcx
	add	%rax, %rbx
	adc	%rdx, %rcx

	mov	24(%rdi), %rax	/* < 2^63-1 */
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^63 */
	add	%rcx, %rbp
	adc	%rdx, %r12	/* CF=0 */


	mov	8(%rdi), %r9
	imul	(%rdi), %r9
	add	%r11, %r9
	neg	%r9

	mov	16(%rsi), %rax
	mul	%r9
	xor	%rcx, %rcx
	add	%r9, %r11
	adc	$0, %r8
	adc	%rax, %rbp
	adc	%rdx, %rcx

	mov	24(%rsi), %rax	/* < 2^63 */
	mul	%r9
	xor	%rbx, %rbx
	add	%rax, %rcx
	adc	%rdx, %rbx	/* < 2^63 */
	xor	%r13, %r13
	add	%rcx, %r12
	adc	%rbx, %r13	/* <= 2^63, so CF=0 */

	mov	8(%rdi), %r9
	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	%r9, %rax
	mul	%r9
	add	%rbx, %r11
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax	/* < 2^63-1 */
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rcx, %r12
	adc	%rdx, %r13	/* CF=0 */


	mov	16(%rdi), %r9
	imul	(%rdi), %r9
	add	%r8, %r9
	neg	%r9

	mov	16(%rsi), %rax
	mul	%r9
	xor	%rbx, %rbx
	xor	%rcx, %rcx
	add	%r9, %r8
	adc	$0, %rbp
	adc	%rax, %rbx
	adc	%rdx, %rcx

	mov	24(%rsi), %rax	/* < 2^63 */
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %r10	/* < 2^63-1 */
	adc	$0, %r10	/* < 2^63 */
	add	%rcx, %r13
	adc	$0, %r10	/* <= 2^63 */


	mov	16(%rdi), %r9
	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	%r9, %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax	/* < 2^63-1 */
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rcx, %r13
	adc	%rdx, %r10	/* CF=0 */


	mov	24(%rdi), %r9
	imul	(%rdi), %r9
	add	%rbp, %r9
	neg	%r9

	mov	16(%rsi), %rax
	mul	%r9
	xor	%rbx, %rbx
	xor	%rcx, %rcx
	add	%r9, %rbp
	adc	$0, %r12
	adc	%rax, %rbx
	adc	%rdx, %rcx

	mov	24(%rsi), %rax	/* < 2^63 */
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	mov	%rdx, %r11	/* < 2^63-1 */
	adc	$0, %r11	/* < 2^63 */
	add	%rcx, %r10
	adc	$0, %r11	/* <= 2^63 */


	mov	24(%rdi), %r9	/* < 2^63-1 */
	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx	/* < 2^63-2 */

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx	/* < 2^63-2 */
	adc	$0, %rbx	/* < 2^63-1 */

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx	/* < 2^63-2 */
	adc	$0, %rcx	/* < 2^63-1 */

	mov	%r9, %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx	/* < 2^62 */
	add	%rcx, %r10
	adc	%rdx, %r11	/* CF=0 */

0:
	mov	%r12, (%rdi)
	mov	%r13, 8(%rdi)
	mov	%r10, 16(%rdi)
	mov	%r11, 24(%rdi)
	sub	$1, %r12
	sbb	$0, %r13
	sbb	16(%rsi), %r10
	sbb	24(%rsi), %r11
	jnc	0b

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
	.globl	_sqrmod256_proth1
	.globl	sqrmod256_proth1
	.align	16

_sqrmod256_proth1:
sqrmod256_proth1:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
#endif
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	24(%rsi), %rsi

	/* %rsi = N[3]
	   %rdi = R  */

	mov	(%rdi), %r9
	mov	%r9, %rax

	mul	%rax
	mov	%rax, %rbx
	mov	%rdx, %r11

	neg	%rax

	mul	%rsi
	mov	%rax, %rbp
	mov	%rdx, %r12

	mov	8(%rdi), %rax
	mul	%r9
	xor	%r8, %r8
	neg	%rbx
	adc	%rax, %r11
	adc	%rdx, %r8

	mov	16(%rdi), %rax
	mul	%r9
	xor	%rcx, %rcx
	add	%rax, %r8
	adc	%rdx, %rcx

	mov	24(%rdi), %rax
	mul	%r9
	add	%rax, %rcx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rcx, %rbp
	adc	%rdx, %r12
	adc	$0, %r13


	mov	8(%rdi), %r9
	imul	(%rdi), %r9
	add	%r11, %r9
	neg	%r9

	mov	%rsi, %rax
	mul	%r9
	xor	%r10, %r10
	add	%r9, %r11
	adc	$0, %r8
	adc	$0, %rbp
	adc	%rax, %r12
	adc	%rdx, %r13
	adc	$0, %r10

	mov	8(%rdi), %r9
	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	%r9, %rax
	mul	%r9
	add	%rbx, %r11
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r8
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r12
	adc	%rdx, %r13
	adc	$0, %r10


	mov	16(%rdi), %r9
	imul	(%rdi), %r9
	add	%r8, %r9
	neg	%r9

	mov	%rsi, %rax
	mul	%r9
	xor	%r11, %r11
	add	%r9, %r8
	adc	$0, %rbp
	adc	$0, %r12
	adc	%rax, %r13
	adc	%rdx, %r10
	adc	$0, %r11

	mov	16(%rdi), %r9
	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %r8
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	%r9, %rax
	mul	%r9
	add	%rcx, %rbp
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	24(%rdi), %rax
	mul	%r9
	add	%rbx, %r12
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r13
	adc	%rdx, %r10
	adc	$0, %r11


	mov	24(%rdi), %r9
	imul	(%rdi), %r9
	add	%rbp, %r9
	neg	%r9

	mov	%rsi, %rax
	mul	%r9

	xor	%r8, %r8	/* carry */
	add	%r9, %rbp
	adc	$0, %r12
	adc	$0, %r13
	adc	%rax, %r10
	adc	%rdx, %r11
	adc	$0, %r8

	mov	24(%rdi), %r9
	mov	(%rdi), %rax
	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r9
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r9
	add	%rcx, %r12
	adc	%rax, %rbx
	mov	%rdx, %rcx
	adc	$0, %rcx

	mov	%r9, %rax
	mul	%r9
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r10
	adc	%rdx, %r11
	adc	$0, %r8

0:
	mov	%r12, (%rdi)
	mov	%r13, 8(%rdi)
	mov	%r10, 16(%rdi)
	mov	%r11, 24(%rdi)
	sub	$1, %r12
	sbb	$0, %r13
	sbb	$0, %r10
	sbb	%rsi, %r11
	sbb	$0, %r8
	jnc	0b

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
