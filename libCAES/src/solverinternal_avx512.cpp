#include "solverinternal.h"
#include "vecmath/vecmath.h"

namespace ECHMET {
namespace CAES {

template <> template <>
void SolverInternal<double, InstructionSet::AVX512>::VectorizedDelogifier<InstructionSet::AVX512>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "AVX512 delogifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
	#if defined(ECHMET_PLATFORM_WIN32) && defined(__x86_64__) && (defined(ECHMET_COMPILER_MINGW) || defined(ECHMET_COMPILER_MSYS))
		m_vecMath.exp10m(src + idx, dst + idx);
	#else
		__m512d _s = M512D(src + idx);

		_s = m_vecMath.exp10m(_s);

		_mm512_store_pd(dst + idx, _s);
	#endif // ECHMET_WIN64_ABI_BUG_CHECK
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.exp10m_single(src[idx]);
}
#ifdef ECHMET_COMPILER_MSVC
template
void SolverInternal<double, InstructionSet::AVX512>::VectorizedDelogifier<InstructionSet::AVX512>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);
#endif // ECHMET_COMPILER_MSVC

template <> template <>
void SolverInternal<double, InstructionSet::AVX512>::VectorizedLogifier<InstructionSet::AVX512>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "AVX512 logifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m512d _s = m_vecMath.mlog10(src + idx);

		_mm512_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.mlog10_single(src[idx]);
}
#ifdef ECHMET_COMPILER_MSVC
template
void SolverInternal<double, InstructionSet::AVX512>::VectorizedLogifier<InstructionSet::AVX512>::operator()(double *ECHMET_RESTRICT_PTR dst, const double *ECHMET_RESTRICT_PTR src);
#endif // ECHMET_COMPILER_MSVC

} // namespace CAES
} // namespace ECHMET
