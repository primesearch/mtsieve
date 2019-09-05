/* mulmod128.S -- (C) Geoffrey Reynolds, March 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void mulmod128(uint64_t *R, const uint64_t *A, const uint64_t *B,
                  const uint64_t *N, uint64_t inv);

   Assign R <-- A*B (mod N) in Montgomery form, where 0 <= A,B < N < 2^128.
   Assumes that A,B are in Montgomery form, N is odd, and N[0]*inv = -1.


   void mulmod128_proth1(uint64_t *R, const uint64_t *A, const uint64_t *B,
                         const uint64_t *N);

   As above, but assumes that N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.text
	.globl	_mulmod128
	.globl	mulmod128
	.p2align 4,,15

_mulmod128:
mulmod128:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
	mov	64(%rsp), %r8
#endif

	mov	(%rdx), %rax
	mov	8(%rdx), %rdx

	mov	(%rcx), %r9
	mov	8(%rcx), %r10

	push	%rbp
	push	%rbx
	push	%r12
	push	%r13

	mov	(%rsi), %rbp
	mov	8(%rsi), %r12

	mov	%rax, (%rdi)
	mov	%rdx, 8(%rdi)

	/* %rdi = B,R
	   %rbp = A[0]
	   %r12 = A[1]
	   %r8  = inv
	   %r9  = N[0]
	   %r10 = N[1] */

	imul	%r8, %rax
	imul	%rbp, %rax
	mov	%rax, %rbx

	mul	%r9
	mov	%rax, %rsi
	mov	%rdx, %r13

	mov	%r10, %rax
	mul	%rbx
	xor	%r11, %r11
	add	%rax, %r13
	adc	%rdx, %r11

	mov	(%rdi), %rax
	mul	%rbp
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%rbp
	add	%rbx, %rsi
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%ebp, %ebp
	add	%rcx, %r13
	adc	%rdx, %r11
	adc	$0, %ebp


	mov	(%rdi), %rax
	mul	%r12
	add	%r13, %rax
	imul	%r8, %rax
	mov	%rax, %rsi

	mul	%r9
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	%r10, %rax
	mul	%rsi
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx
	add	%rcx, %r11
	adc	%rdx, %rbp	/* CF=0 */

	mov	(%rdi), %rax
	mul	%r12
	mov	%rax, %rbx
	mov	%rdx, %rcx

	mov	8(%rdi), %rax
	mul	%r12
	add	%rbx, %r13
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%esi, %esi
	add	%r11, %rcx
	adc	%rbp, %rdx
	adc	$0, %esi

0:	mov	%rcx, (%rdi)
	mov	%rdx, 8(%rdi)
	sub	%r9, %rcx
	sbb	%r10, %rdx
	sbb	$0, %esi
	jnc	0b

	pop	%r13
	pop	%r12
	pop	%rbx
	pop	%rbp
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	ret


	.text
	.globl	_mulmod128_proth1
	.globl	mulmod128_proth1
	.p2align 4,,15

_mulmod128_proth1:
mulmod128_proth1:
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi
	mov	%rdx, %rsi
	mov	%r8, %rdx
	mov	%r9, %rcx
#endif
	mov	(%rdx), %rax
	mov	8(%rdx), %rdx
	mov	8(%rcx), %r10

	push	%rbp
	push	%rbx

	mov	(%rsi), %rbp
	mov	8(%rsi), %r8

	mov	%rax, (%rdi)
	mov	%rdx, 8(%rdi)

	/* %rdi = B,R
	   %rbp = A[0]
	   %r8  = A[1]
	   %r10 = N[1] */

	mul	%rbp
	neg	%rax
	mov	%rdx, %rcx
	mov	%rax, %rbx

	mul	%r10
	mov	%rax, %r9
	mov	%rdx, %r11

	mov	8(%rdi), %rax
	mul	%rbp
	neg	%rbx		/* set carry if nonzero */
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%ebp, %ebp
	add	%rcx, %r9
	adc	%rdx, %r11
	adc	$0, %ebp


	mov	(%rdi), %rax
	mul	%r8
	mov	%rax, %rbx
	mov	%rdx, %rcx
	add	%r9, %rax
	neg	%rax
	mov	%rax, %rsi

	mul	%r10
	add	%rsi, %r9
	adc	%rax, %r11
	adc	%rdx, %rbp	/* CF=0 */

	mov	8(%rdi), %rax
	mul	%r8
	add	%rbx, %r9
	adc	%rax, %rcx
	adc	$0, %rdx
	xor	%esi, %esi
	add	%r11, %rcx
	adc	%rbp, %rdx
	adc	$0, %esi

0:	mov	%rcx, (%rdi)
	mov	%rdx, 8(%rdi)
	sub	$1, %rcx
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
