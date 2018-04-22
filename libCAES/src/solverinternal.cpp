#include "solverinternal.h"
#include "vecmath/vecmath.h"

namespace ECHMET {
namespace CAES {

template <> template <>
__attribute__((target("sse2")))
void SolverInternal<double, InstructionSet::SSE2>::VectorizedDelogifier<InstructionSet::SSE2>::operator()(double  *__restrict__ dst, const double *__restrict__ src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "SSE2 delogifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m128d _s = *(__m128d *)(src + idx);

		_s = m_vecMath.exp10m(_s);

		_mm_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.exp10m_single(src[idx]);
}

template <> template <>
__attribute__((target("sse2")))
void SolverInternal<double, InstructionSet::SSE2>::VectorizedLogifier<InstructionSet::SSE2>::operator()(double  *__restrict__ dst, const double *__restrict__ src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "SSE2 logifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m128d _s = m_vecMath.mlog10(src + idx);

		_mm_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.mlog10_single(src[idx]);
}

template <> template <>
__attribute__((target("avx")))
void SolverInternal<double, InstructionSet::AVX>::VectorizedDelogifier<InstructionSet::AVX>::operator()(double  *__restrict__ dst, const double *__restrict__ src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "AVX delogifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m256d _s = *(__m256d *)(src + idx);

		_s = m_vecMath.exp10m(_s);

		_mm256_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.exp10m_single(src[idx]);
}

template <> template <>
__attribute__((target("avx")))
void SolverInternal<double, InstructionSet::AVX>::VectorizedLogifier<InstructionSet::AVX>::operator()(double  *__restrict__ dst, const double *__restrict__ src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "AVX logifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m256d _s = m_vecMath.mlog10(src + idx);

		_mm256_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.mlog10_single(src[idx]);
}

template <> template <>
__attribute__((target("fma")))
void SolverInternal<double, InstructionSet::FMA3>::VectorizedDelogifier<InstructionSet::FMA3>::operator()(double  *__restrict__ dst, const double *__restrict__ src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "FMA3 delogifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m256d _s = *(__m256d *)(src + idx);

		_s = m_vecMath.exp10m(_s);

		_mm256_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.exp10m_single(src[idx]);
}

template <> template <>
__attribute__((target("fma")))
void SolverInternal<double, InstructionSet::FMA3>::VectorizedLogifier<InstructionSet::FMA3>::operator()(double  *__restrict__ dst, const double *__restrict__ src)
{
	ECHMET_DEBUG_CODE(fprintf(stderr, "FMA3 logifier, NBlock: %zu, total: %zu\n", NBlock, N));

	size_t idx = 0;
	for (; idx < NBlock; idx += blockSize) {
		__m256d _s = m_vecMath.mlog10(src + idx);

		_mm256_store_pd(dst + idx, _s);
	}

	for (; idx < N; idx++)
		dst[idx] = m_vecMath.mlog10_single(src[idx]);
}

} // namespace CAES
} // namespace ECHMET
