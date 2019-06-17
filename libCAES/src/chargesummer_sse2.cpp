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
		__m128d chg = M128D(m_charges + idx);

		zVec = _mm_add_pd(zVec, _mm_mul_pd(conc, chg));
	}

	_mm_store_pd(vz, zVec);
	z = vz[0] + vz[1];

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

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
		__m128d chg = M128D(m_charges + idx);

		zVec = _mm_add_pd(zVec, _mm_mul_pd(conc, chg));
	}

	_mm_store_pd(vz, zVec);
	z = vz[0] + vz[1];

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

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
		__m128d chg = M128D(m_charges + idx);

		zVec = _mm_add_pd(zVec, _mm_mul_pd(conc, chg));
		dZVec = _mm_add_pd(dZVec, _mm_mul_pd(dConc, chg));
	}

	_mm_store_pd(vz, zVec);
	_mm_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1];
	dZ = vdZ[0] + vdZ[1];

	for (; idx < m_N; idx++) {
		z += m_charges[idx] * icConcs[idx];
		dZ += m_charges[idx] * dIcConcsdH[idx];
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
		__m128d chg = M128D(m_charges + idx);

		zVec = _mm_add_pd(zVec, _mm_mul_pd(conc, chg));
		dZVec = _mm_add_pd(dZVec, _mm_mul_pd(dConc, chg));
	}

	_mm_store_pd(vz, zVec);
	_mm_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1];
	dZ = vdZ[0] + vdZ[1];

	for (; idx < m_N; idx++) {
		z += m_charges[idx] * icConcs[idx];
		dZ += m_charges[idx] * dIcConcsdH[idx];
	}
}

template <>
double ChargeSummer<double, InstructionSet::SSE2, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	VD vIs;
	__m128d isVec = M128D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m128d conc = M128D(icConcs + idx);
		__m128d chgSq = M128D(m_chargesSquared + idx);

		isVec = _mm_add_pd(isVec, _mm_mul_pd(conc, chgSq));
	}

	_mm_store_pd(vIs, isVec);

	is = vIs[0] + vIs[1];

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}

template <>
double ChargeSummer<double, InstructionSet::SSE2, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	VD vIs;
	__m128d isVec = M128D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m128d conc = M128D(icConcs + idx);
		__m128d chgSq = M128D(m_chargesSquared + idx);

		isVec = _mm_add_pd(isVec, _mm_mul_pd(conc, chgSq));
	}

	_mm_store_pd(vIs, isVec);

	is = vIs[0] + vIs[1];

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}

} // namespace CAES
} // namespace ECHMET
