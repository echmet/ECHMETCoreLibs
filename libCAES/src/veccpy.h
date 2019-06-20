#ifndef ECHMET_CAES_VECCPY_H
#define ECHMET_CAES_VECCPY_H

#include <cstddef>

#define ECHMET_IMPORT_INTERNAL
#include <echmetmodule.h>
#undef ECHMET_IMPORT_INTERNAL

namespace ECHMET {
namespace CAES {

#ifdef ECHMET_COMPILER_GCC_LIKE

inline
void smlcpy(double *const dst, const double *const src, const size_t len)
{
	asm volatile(
		"movq $0, %%rcx;"
		"smlcpy_loop%=: ;"
		"movq (%1,%%rcx,8), %%rax;"
		"movq %%rax, (%0,%%rcx,8);"
		"add $1, %%rcx;"
		"cmpq %2, %%rcx;"
		"jl smlcpy_loop%=;"
		:
		: "r"(dst), "r"(src), "r"(len)
		: "%rcx", "%rax" );
}

inline
void midcpy(double *const dst, const double *const src, const size_t len)
{
	const size_t blockLen{len- (len % 2)};

	asm volatile(
		"movq $0, %%rcx;"
		"midcpy_loop%=: ;"
		"movq (%1,%%rcx,8), %%rax;"
		"movq 8(%1,%%rcx,8), %%rbx;"
		"movq %%rax, (%0,%%rcx,8);"
		"movq %%rbx, 8(%0,%%rcx,8);"
		"add $2, %%rcx;"
		"cmpq %2, %%rcx;"
		"jl midcpy_loop%=;"
		"cmpq %3, %%rcx;"
		"je midcpy_out%=; "
		"movq (%1,%%rcx,8), %%rax;"
		"movq %%rax, (%0,%%rcx,8);"
		"midcpy_out%=: ;"
		:
		: "r"(dst), "r"(src), "r"(blockLen), "r"(len)
		: "%rcx", "%rax", "%rbx" );
}

inline
void dbl_vec_cpy(double *const ECHMET_RESTRICT_PTR dst, const double *const ECHMET_RESTRICT_PTR src, const size_t len)
{
	if (len < 3)
		smlcpy(dst, src, len);
	else if (len < 32)
		midcpy(dst, src, len);
	else
		memcpy(dst, src, len * 8);
}

#else

inline
void dbl_vec_cpy(double *const ECHMET_RESTRICT_PTR dst, const double *const ECHMET_RESTRICT_PTR src, const size_t len)
{
	memcpy(dst, src, len * sizeof(double);
}

#endif // ECHMET_COMPILER_GCC_LIKE

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_VECCPY_H
