#include "chargesummer.h"

namespace ECHMET {
namespace CAES {

using VD = typename VecMath<InstructionSet::AVX>::VD;

static const VD ZERO = {0, 0, 0, 0};

template <>
double ChargeSummer<double, InstructionSet::AVX, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	VD vz;
	__m256d zVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chg = M256D(m_chargesArray + idx);

		conc = _mm256_mul_pd(conc, chg);
		zVec = _mm256_add_pd(zVec, conc);
	}

	_mm256_store_pd(vz, zVec);
	z = vz[0] + vz[1];

	for (; idx < m_N; idx++)
		z += m_chargesArray[idx] * icConcs[idx];

	return z;
}

template <>
double ChargeSummer<double, InstructionSet::AVX, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	VD vz;
	__m256d zVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chg = M256D(m_chargesArray + idx);

		conc = _mm256_mul_pd(conc, chg);
		zVec = _mm256_add_pd(zVec, conc);
	}

	_mm256_store_pd(vz, zVec);
	z = vz[0] + vz[1];

	for (; idx < m_N; idx++)
		z += m_chargesArray[idx] * icConcs[idx];

	return z;
}

template <>
void ChargeSummer<double, InstructionSet::AVX, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								   const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								   double &z, double &dZ) noexcept
{
	VD vz;
	VD vdZ;
	__m256d zVec = M256D(ZERO);
	__m256d dZVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d dConc = M256D(dIcConcsdH + idx);
		__m256d chg = M256D(m_chargesArray + idx);

		conc = _mm256_mul_pd(conc, chg);
		dConc = _mm256_mul_pd(dConc, chg);

		zVec = _mm256_add_pd(zVec, conc);
		dZVec = _mm256_add_pd(dZVec, dConc);
	}

	_mm256_store_pd(vz, zVec);
	_mm256_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1];
	dZ = vdZ[0] + vdZ[1];

	for (; idx < m_N; idx++) {
		z += m_chargesArray[idx] * icConcs[idx];
		dZ += m_chargesArray[idx] * dIcConcsdH[idx];
	}
}

template <>
void ChargeSummer<double, InstructionSet::AVX, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept
{
	VD vz;
	VD vdZ;
	__m256d zVec = M256D(ZERO);
	__m256d dZVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d dConc = M256D(dIcConcsdH + idx);
		__m256d chg = M256D(m_chargesArray + idx);

		conc = _mm256_mul_pd(conc, chg);
		dConc = _mm256_mul_pd(dConc, chg);

		zVec = _mm256_add_pd(zVec, conc);
		dZVec = _mm256_add_pd(dZVec, dConc);
	}

	_mm256_store_pd(vz, zVec);
	_mm256_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1];
	dZ = vdZ[0] + vdZ[1];

	for (; idx < m_N; idx++) {
		z += m_chargesArray[idx] * icConcs[idx];
		dZ += m_chargesArray[idx] * dIcConcsdH[idx];
	}
}

} // namespace CAES
} // namespace ECHMET
