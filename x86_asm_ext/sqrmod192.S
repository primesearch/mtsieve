/* sqrmod192.S -- (C) Geoffrey Reynolds, March 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void sqrmod192(uint64_t *R, const uint64_t *N, uint64_t inv);

   Assign R <-- R^2 (mod N) in Montgomery form, where 0 <= R < N < 2^192.
   Assumes that R is in Montgomery form, N is odd, and N[0]*inv = -1.


   void sqrmod192_proth0(uint64_t *R, const uint64_t *N);

   As sqrmod192, but assumes that N[2] < 2^63, N[0] = 1.


   void sqrmod192_proth1(uint64_t *R, const uint64_t *N);

   As sqrmod192, but assumes that N[1] = 0, N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.text
	.globl	_sqrmod192
	.globl	sqrmod192
	.p2align 4,,15

_sqrmod192:
sqrmod192:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
#endif
	mov	%rdx, %r8
	mov	%rsi, %r9

	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	/* %r8  = inv
	   %r9  = N
	   %rdi = R */

	mov	(%rdi), %rax
	mul	%rax
	mul	%r8
	mov	%rax, %r10

	mov	(%r9), %rax
	mul	%r10
	mov	%rax, %r12
	mov	%rdx, %r14

	mov	8(%r9), %rax
	mul	%r10
	add	%rax, %r14
	mov	%rdx, %rbp
	adc	$0, %rbp

	mov	16(%r9), %rax
	mul	%r10
	add	%rax, %rbp
	mov	%rdx, %r11
	adc	$0, %r11

	mov	(%rdi), %rax
	mul	%rax
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mulq	(%rdi)
	add	%rbx, %r12
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mulq	(%rdi)
	add	%rcx, %r14
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r12, %r12
	add	%rbx, %rbp
	adc	%rdx, %r11
	adc	$0, %r12


	mov	8(%rdi), %rax
	mulq	(%rdi)
	add	%r14, %rax
	mul	%r8
	mov	%rax, %r10

	mov	(%r9), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r9), %rax
	mul	%r10
	add	%rbx, %r14
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r9), %rax
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r13, %r13
	add	%rbx, %r11
	adc	%rdx, %r12
	adc	$0, %r13

	mov	8(%rdi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %r14
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r10
	add	%rcx, %rbp
	adc	%rax, %rbx
	adc	$0, %rdx
	add	%rbx, %r11
	adc	%rdx, %r12
	adc	$0, %r13


	mov	16(%rdi), %rax
	mulq	(%rdi)
	add	%rbp, %rax
	mul	%r8
	mov	%rax, %r10

	mov	(%r9), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%r9), %rax
	mul	%r10
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%r9), %rax
	mul	%r10
	add	%rcx, %r11
	adc	%rax, %rbx
	adc	$0, %rdx
	xor	%r8, %r8
	add	%rbx, %r12
	adc	%rdx, %r13
	adc	$0, %r8

	mov	16(%rdi), %r10
	mov	(%rdi), %rax
	mul	%r10
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r10
	add	%rbx, %rbp
	adc	%rax, %rcx
	mov	%rdx, %rbx
	adc	$0, %rbx

	mov	16(%rdi), %rax
	mul	%r10
	add	%r11, %rcx
	adc	%rbx, %rax
	adc	$0, %rdx
	add	%r12, %rax
	adc	%r13, %rdx
	adc	$0, %r8

0:	mov	%rcx, (%rdi)
	mov	%rax, 8(%rdi)
	mov	%rdx, 16(%rdi)
	sub	(%r9), %rcx
	sbb	8(%r9), %rax
	sbb	16(%r9), %rdx
	sbb	$0, %r8
	jnc	0b

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
	.globl	_sqrmod192_proth0
	.globl	sqrmod192_proth0
	.p2align 4,,15

_sqrmod192_proth0:
sqrmod192_proth0:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
#endif
	mov	8(%rsi), %r9
	mov	16(%rsi), %r8

	push	%rbx
	push	%rbp

	/* %r9  = N[1]
	   %r8  = N[2] < 2^63
	   %rdi = R */

	mov	(%rdi), %rax
	xor	%ecx, %ecx
	mul	%rax
	neg	%rax
	mov	%rax, %rbx
	adc	%rdx, %rcx

	mul	%r9
	xor	%ebp, %ebp
	add	%rax, %rcx
	adc	%rdx, %rbp

	mov	%r8, %rax	/* < 2^63 */
	mul	%rbx
	mov	(%rdi), %r10
	xor	%r11, %r11
	add	%rax, %rbp
	adc	%rdx, %r11	/* < 2^63 */

	mov	8(%rdi), %rax
	mul	%r10
	mov	%rax, %rsi
	add	%rax, %rcx
	adc	%rdx, %rbp
	adc	$0, %r11	/* <= 2^63 */

	mov	16(%rdi), %rax	/* < 2^63-1 */
	mul	%r10
	add	%rax, %rbp
	adc	%rdx, %r11	/* CF=0 */

	lea	(%rcx,%rsi), %rax
	neg	%rax
	mov	%rax, %r10

	mul	%r9
	xor	%ebx, %ebx
	add	%r10, %rcx
	adc	%rax, %rbp
	adc	%rdx, %r11
	adc	$0, %rbx	/* < 2 */

	mov	%r8, %rax	/* < 2^63 */
	mul	%r10
	mov	8(%rdi), %r10
	add	%rax, %r11
	adc	%rdx, %rbx	/* < 2^63 */

	mov	(%rdi), %rax
	xor	%esi, %esi
	mul	%r10
	add	%rax, %rcx
	adc	%rdx, %rsi	/* CF=0 */

	mov	%r10, %rax
	mul	%rax
	add	%rsi, %rbp
	adc	$0, %rdx	/* CF=0 */
	add	%rax, %rbp
	adc	%rdx, %r11
	adc	$0, %rbx	/* <= 2^63 */

	mov	16(%rdi), %rcx	/* < 2^63-1 */
	mov	%rcx, %rax
	mul	%r10
	add	%rax, %r11
	adc	%rdx, %rbx	/* CF=0 */

	mov	(%rdi), %rax
	imul	%rcx, %rax	/* %rsi < 2^63-1 */
	add	%rbp, %rax
	neg	%rax
	mov	%rax, %r10

	mul	%r9
	xor	%ecx, %ecx
	add	%r10, %rbp
	adc	%rax, %r11
	adc	%rdx, %rbx
	adc	$0, %rcx	/* < 2 */

	mov	%r8, %rax	/* < 2^63 */
	mul	%r10
	mov	16(%rdi), %r10	/* < 2^63-1 */
	add	%rax, %rbx
	adc	%rdx, %rcx	/* <= 2^63 */

	mov	(%rdi), %rax
	mul	%r10
	xor	%esi, %esi
	add	%rax, %rbp
	adc	%rdx, %rsi	/* < 2^63-1 */

	mov	8(%rdi), %rax
	mul	%r10		/* < 2^63-1 */
	add	%r11, %rsi
	adc	$0, %rdx	/* < 2^63-1 */
	add	%rax, %rsi
	adc	%rdx, %rbx
	adc	$0, %rcx	/* <= 2^63+1 */

	mov	%r10, %rax	/* < 2^63-1 */
	mul	%rax
	add	%rax, %rbx
	adc	%rdx, %rcx	/* CF=0 */

0:	mov	%rsi, (%rdi)
	mov	%rbx, 8(%rdi)
	mov	%rcx, 16(%rdi)
	sub	$1, %rsi
	sbb	%r9, %rbx
	sbb	%r8, %rcx
	jnc	0b

	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret




	.text
	.globl	_sqrmod192_proth1
	.globl	sqrmod192_proth1
	.p2align 4,,15

_sqrmod192_proth1:
sqrmod192_proth1:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
#endif
	mov	(%rdi), %rax
	mov	8(%rdi), %r10
	mov	16(%rdi), %r8
	mov	16(%rsi), %r9

	push	%rbx
	push	%rbp
	push	%r12

	/* %r9  = N[2]
	   %r10 = R[1]
	   %r8  = R[2]
	   %rdi = R  */

	mul	%rax
	xor	%esi, %esi
	neg	%rax
	adc	%rdx, %rsi	/* CF=0 */

	mul	%r9
	mov	%rax, %rbp
	mov	%rdx, %r11	/* < 2^64-1 */

	mov	(%rdi), %rax
	mul	%r10
	xor	%ecx, %ecx
	add	%rax, %rsi
	adc	%rdx, %rcx	/* CF=0 */
	mov	%rax, %rbx

	mov	(%rdi), %rax
	xor	%r12, %r12
	mul	%r8
	add	%rcx, %rbp
	adc	$0, %r11	/* CF=0 */
	add	%rax, %rbp
	adc	%rdx, %r11
	adc	$0, %r12	/* < 2 */

	lea	(%rbx,%rsi), %rax
	neg	%rax
	mov	%rax, %rbx

	mul	%r9
	add	%rbx, %rsi
	adc	$0, %rbp
	adc	%rax, %r11
	adc	%rdx, %r12	/* CF=0 */

	mov	(%rdi), %rax
	mul	%r10
	xor	%ecx, %ecx
	add	%rax, %rsi
	adc	%rdx, %rcx

	mov	%r10, %rax
	mul	%rax
	xor	%ebx, %ebx
	add	%rax, %rcx
	adc	%rdx, %rbx	/* CF=0 */

	mov	%r10, %rax
	mul	%r8
	add	%rcx, %rbp
	adc	%rbx, %r11
	adc	$0, %rdx	/* CF=0 */
	xor	%esi, %esi
	add	%rax, %r11
	adc	%rdx, %r12
	adc	$0, %esi	/* < 2 */

	mov	(%rdi), %rax
	imul	%r8, %rax
	add	%rbp, %rax
	neg	%rax

	xor	%ecx, %ecx
	add	%rax, %rbp
	adc	$0, %rcx	/* < 2 */

	mul	%r9
	xor	%r10, %r10
	add	%rcx, %r11
	adc	%rax, %r12
	adc	%rdx, %rsi
	adc	$0, %r10	/* < 2 */

	mov	(%rdi), %rax
	mul	%r8
	xor	%ecx, %ecx
	add	%rax, %rbp
	adc	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r8
	xor	%ebx, %ebx
	add	%rax, %rcx
	adc	%rdx, %rbx	/* CF=0 */

	mov	%r8, %rax
	mul	%rax
	add	%r11, %rcx
	adc	%r12, %rbx
	adc	$0, %rdx	/* CF=0 */
	add	%rax, %rbx
	adc	%rdx, %rsi
	adc	$0, %r10

0:	mov	%rcx, (%rdi)
	mov	%rbx, 8(%rdi)
	mov	%rsi, 16(%rdi)
	sub	$1, %rcx
	sbb	$0, %rbx
	sbb	%r9, %rsi
	sbb	$0, %r10
	jnc	0b

	pop	%r12
	pop	%rbp
	pop	%rbx
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret
