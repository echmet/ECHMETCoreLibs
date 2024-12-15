#include "chargesummer.h"

namespace ECHMET {
namespace CAES {

using VD = typename VecMath<InstructionSet::AVX512>::VD;

template <>
double ChargeSummer<double, InstructionSet::AVX512, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	__m512d zVec = _mm512_setzero_pd();

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m512d conc = M512D(icConcs + idx);
		__m512d chg = M512D(m_charges + idx);

		zVec = _mm512_fmadd_pd(conc, chg, zVec);
	}

	z = _mm512_reduce_add_pd(zVec);

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

	return z;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::AVX512, false>::calc(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
double ChargeSummer<double, InstructionSet::AVX512, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	__m512d zVec = _mm512_setzero_pd();

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m512d conc = M512D(icConcs + idx);
		__m512d chg = M512D(m_charges + idx);

		zVec = _mm512_fmadd_pd(conc, chg, zVec);
	}

	z = _mm512_reduce_add_pd(zVec);

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

	return z;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::AVX512, true>::calc(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
void ChargeSummer<double, InstructionSet::AVX512, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								   const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								   double &z, double &dZ) noexcept
{
	__m512d zVec = _mm512_setzero_pd();
	__m512d dZVec = _mm512_setzero_pd();

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m512d conc = M512D(icConcs + idx);
		__m512d dConc = M512D(dIcConcsdH + idx);
		__m512d chg = M512D(m_charges + idx);

		zVec = _mm512_fmadd_pd(conc, chg, zVec);
		dZVec = _mm512_fmadd_pd(dConc, chg, dZVec);
	}

	z = _mm512_reduce_add_pd(zVec);
	dZ = _mm512_reduce_add_pd(dZVec);

	for (; idx < m_N; idx++) {
		z += m_charges[idx] * icConcs[idx];
		dZ += m_charges[idx] * dIcConcsdH[idx];
	}
}
#ifdef ECHMET_COMPILER_MSVC
template
void ChargeSummer<double, InstructionSet::AVX512, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR,
	const double *const ECHMET_RESTRICT_PTR,
	double &, double &) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
void ChargeSummer<double, InstructionSet::AVX512, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept
{
	__m512d zVec = _mm512_setzero_pd();
	__m512d dZVec = _mm512_setzero_pd();

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m512d conc = M512D(icConcs + idx);
		__m512d dConc = M512D(dIcConcsdH + idx);
		__m512d chg = M512D(m_charges + idx);

		zVec = _mm512_fmadd_pd(conc, chg, zVec);
		dZVec = _mm512_fmadd_pd(dConc, chg, dZVec);
	}

	z = _mm512_reduce_add_pd(zVec);
	dZ = _mm512_reduce_add_pd(dZVec);

	for (; idx < m_N; idx++) {
		z += m_charges[idx] * icConcs[idx];
		dZ += m_charges[idx] * dIcConcsdH[idx];
	}
}
#ifdef ECHMET_COMPILER_MSVC
template
void ChargeSummer<double, InstructionSet::AVX512, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR,
	const double *const ECHMET_RESTRICT_PTR,
	double &, double &) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
double ChargeSummer<double, InstructionSet::AVX512, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	__m512d isVec = _mm512_setzero_pd();

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m512d conc = M512D(icConcs + idx);
		__m512d chgSq = M512D(m_chargesSquared + idx);

		isVec = _mm512_fmadd_pd(conc, chgSq, isVec);
	}

	is = _mm512_reduce_add_pd(isVec);

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::AVX512, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
double ChargeSummer<double, InstructionSet::AVX512, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	__m512d isVec = _mm512_setzero_pd();

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m512d conc = M512D(icConcs + idx);
		__m512d chgSq = M512D(m_chargesSquared + idx);

		isVec = _mm512_fmadd_pd(conc, chgSq, isVec);
	}

	is = _mm512_reduce_add_pd(isVec);

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::AVX512, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

} // namespace CAES
} // namespace ECHMET
