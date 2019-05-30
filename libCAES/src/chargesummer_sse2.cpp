#include "chargesummer.h"

namespace ECHMET {
namespace CAES {

using VD = typename VecMath<InstructionSet::SSE2>::VD;

static const VD ZERO = {0, 0};

template <>
double ChargeSummer<double, InstructionSet::SSE2, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	VD vz;
	__m128d zVec = M128D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m128d conc = M128D(icConcs + idx);
		__m128d chg = M128D(m_chargesArray + idx);

		conc = _mm_mul_pd(conc, chg);
		zVec = _mm_add_pd(zVec, conc);
	}

	_mm_store_pd(vz, zVec);
	z = vz[0] + vz[1];

	for (; idx < m_N; idx++)
		z += m_chargesArray[idx] * icConcs[idx];

	return z;
}

template <>
double ChargeSummer<double, InstructionSet::SSE2, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	VD vz;
	__m128d zVec = M128D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m128d conc = M128D(icConcs + idx);
		__m128d chg = M128D(m_chargesArray + idx);

		conc = _mm_mul_pd(conc, chg);
		zVec = _mm_add_pd(zVec, conc);
	}

	_mm_store_pd(vz, zVec);
	z = vz[0] + vz[1];

	for (; idx < m_N; idx++)
		z += m_chargesArray[idx] * icConcs[idx];

	return z;
}

template <>
void ChargeSummer<double, InstructionSet::SSE2, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								   const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								   double &z, double &dZ) noexcept
{
	VD vz;
	VD vdZ;
	__m128d zVec = M128D(ZERO);
	__m128d dZVec = M128D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m128d conc = M128D(icConcs + idx);
		__m128d dConc = M128D(dIcConcsdH + idx);
		__m128d chg = M128D(m_chargesArray + idx);

		conc = _mm_mul_pd(conc, chg);
		dConc = _mm_mul_pd(dConc, chg);

		zVec = _mm_add_pd(zVec, conc);
		dZVec = _mm_add_pd(dZVec, dConc);
	}

	_mm_store_pd(vz, zVec);
	_mm_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1];
	dZ = vdZ[0] + vdZ[1];

	for (; idx < m_N; idx++) {
		z += m_chargesArray[idx] * icConcs[idx];
		dZ += m_chargesArray[idx] * dIcConcsdH[idx];
	}
}

template <>
void ChargeSummer<double, InstructionSet::SSE2, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept
{
	VD vz;
	VD vdZ;
	__m128d zVec = M128D(ZERO);
	__m128d dZVec = M128D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m128d conc = M128D(icConcs + idx);
		__m128d dConc = M128D(dIcConcsdH + idx);
		__m128d chg = M128D(m_chargesArray + idx);

		conc = _mm_mul_pd(conc, chg);
		dConc = _mm_mul_pd(dConc, chg);

		zVec = _mm_add_pd(zVec, conc);
		dZVec = _mm_add_pd(dZVec, dConc);
	}

	_mm_store_pd(vz, zVec);
	_mm_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1];
	dZ = vdZ[0] + vdZ[1];

	for (; idx < m_N; idx++) {
		z += m_chargesArray[idx] * icConcs[idx];
		dZ += m_chargesArray[idx] * dIcConcsdH[idx];
	}
}

} // namespace CAES
} // namespace ECHMET