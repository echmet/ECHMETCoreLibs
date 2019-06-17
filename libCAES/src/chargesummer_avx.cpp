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
		__m256d chg = M256D(m_charges + idx);

		zVec = _mm256_add_pd(zVec, _mm256_mul_pd(conc, chg));
	}

	_mm256_store_pd(vz, zVec);
	z = vz[0] + vz[1] + vz[2] + vz[3];

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

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
		__m256d chg = M256D(m_charges + idx);

		zVec = _mm256_add_pd(zVec, _mm256_mul_pd(conc, chg));
	}

	_mm256_store_pd(vz, zVec);
	z = vz[0] + vz[1] + vz[2] + vz[3];

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

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
		__m256d chg = M256D(m_charges + idx);

		zVec = _mm256_add_pd(zVec, _mm256_mul_pd(conc, chg));
		dZVec = _mm256_add_pd(dZVec, _mm256_mul_pd(dConc, chg));
	}

	_mm256_store_pd(vz, zVec);
	_mm256_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1] + vz[2] + vz[3];
	dZ = vdZ[0] + vdZ[1] + vdZ[2] + vdZ[3];

	for (; idx < m_N; idx++) {
		z += m_charges[idx] * icConcs[idx];
		dZ += m_charges[idx] * dIcConcsdH[idx];
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
		__m256d chg = M256D(m_charges + idx);

		zVec = _mm256_add_pd(zVec, _mm256_mul_pd(conc, chg));
		dZVec = _mm256_add_pd(dZVec, _mm256_mul_pd(dConc, chg));
	}

	_mm256_store_pd(vz, zVec);
	_mm256_store_pd(vdZ, dZVec);

	z = vz[0] + vz[1] + vz[2] + vz[3];
	dZ = vdZ[0] + vdZ[1] + vdZ[2] + vdZ[3];

	for (; idx < m_N; idx++) {
		z += m_charges[idx] * icConcs[idx];
		dZ += m_charges[idx] * dIcConcsdH[idx];
	}
}

template <>
double ChargeSummer<double, InstructionSet::AVX, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	VD vIs;
	__m256d isVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chgSq = M256D(m_chargesSquared + idx);

		isVec = _mm256_add_pd(isVec, _mm256_mul_pd(conc, chgSq));
	}

	_mm256_store_pd(vIs, isVec);

	is = vIs[0] + vIs[1] + vIs[2] + vIs[3];

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}

template <>
double ChargeSummer<double, InstructionSet::AVX, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	VD vIs;
	__m256d isVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chgSq = M256D(m_chargesSquared + idx);

		isVec = _mm256_add_pd(isVec, _mm256_mul_pd(conc, chgSq));
	}

	_mm256_store_pd(vIs, isVec);

	is = vIs[0] + vIs[1] + vIs[2] + vIs[3];

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}

} // namespace CAES
} // namespace ECHMET
