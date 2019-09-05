/* redc.S -- (C) Geoffrey Reynolds, December 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void redc_proth0(uint64_t *X, const uint64_t *N, unsigned int n);

   Assumes N = k*2^n+1, k < 2^64, n > 192.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.text
	.globl	_redc_proth0
	.globl	redc_proth0
	.align	16

_redc_proth0:
redc_proth0:
	push	%rbx
#ifdef _WIN64
	push	%rdi
	push	%rsi
	mov	%rcx, %rdi		/* X */
	lea	-16(%rdx,%r8,8), %rsi	/* N+nn-2 */
	mov	%r8, %r11		/* outer counter nn */
	lea	-3(%r8), %r10		/* inner counter nn-3 */
#else
	lea	-16(%rsi,%rdx,8), %rsi	/* N+nn-2 */
	mov	%rdx, %r11		/* outer counter nn */
	lea	-3(%rdx), %r10		/* inner counter nn-3 */
#endif

	.align	16
outer:
	mov	(%rdi), %r9
	mov	%rdi, %r8
	mov	%r10d, %ebx
	xor	%ecx, %ecx		/* initial carry */
	movq	$0, (%rdi)

	/* j = 0 */
	neg	%r9
	setc	%cl
	lea	8(%rdi), %rdi

	test	$1, %ebx
	jz	inner

	/* j = 1 */
	add	%rcx, (%rdi)
	setc	%cl
	lea	8(%rdi), %rdi
	dec	%ebx

	/* j < nn-2 */

	.align	32
inner:
	mov	(%rdi), %rax
	mov	8(%rdi), %rdx
	add	%rcx, %rax
	adc	$0, %rdx
	mov	%rax, (%rdi)
	mov	%rdx, 8(%rdi)
	setc	%cl

	lea	16(%rdi), %rdi
	sub	$2, %ebx
	ja	inner

	/* j = nn-2 */

	movq	(%rsi), %rax		/* N[nn-2] */
	mulq	%r9
	addq	%rcx, %rax
	adcq	$0, %rdx
	addq	%rax, (%rdi)
	adcq	$0, %rdx
	movq	%rdx, %rcx

	/* j = nn-1 */

	movq	8(%rsi), %rax		/* N[nn-1] */
	mulq	%r9
	addq	%rcx, %rax
	adcq	$0, %rdx
	addq	%rax, 8(%rdi)
	adcq	$0, %rdx
	movq	%rdx, %rcx


	mov	%r8, %rdi
	mov	%rcx, (%rdi)		/* save final carry */

	lea	8(%rdi), %rdi
	sub	$1, %r11
	ja	outer

#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	pop	%rbx
	ret
