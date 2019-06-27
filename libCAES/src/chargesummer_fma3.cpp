#include "chargesummer.h"

namespace ECHMET {
namespace CAES {

using VD = typename VecMath<InstructionSet::FMA3>::VD;

static const VD ZERO = {0, 0, 0, 0};

template <>
double ChargeSummer<double, InstructionSet::FMA3, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	VD vz;
	__m256d zVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chg = M256D(m_charges + idx);

		zVec = _mm256_fmadd_pd(conc, chg, zVec);
	}

	_mm256_store_pd(vz, zVec);
	z = vz[0] + vz[1] + vz[2] + vz[3];

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

	return z;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::FMA3, false>::calc(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
double ChargeSummer<double, InstructionSet::FMA3, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double z;
	VD vz;
	__m256d zVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chg = M256D(m_charges + idx);

		zVec = _mm256_fmadd_pd(conc, chg, zVec);
	}

	_mm256_store_pd(vz, zVec);
	z = vz[0] + vz[1] + vz[2] + vz[3];

	for (; idx < m_N; idx++)
		z += m_charges[idx] * icConcs[idx];

	return z;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::FMA3, true>::calc(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
void ChargeSummer<double, InstructionSet::FMA3, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
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

		zVec = _mm256_fmadd_pd(conc, chg, zVec);
		dZVec = _mm256_fmadd_pd(dConc, chg, dZVec);
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
#ifdef ECHMET_COMPILER_MSVC
template
void ChargeSummer<double, InstructionSet::FMA3, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR,
	const double *const ECHMET_RESTRICT_PTR,
	double &, double &) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
void ChargeSummer<double, InstructionSet::FMA3, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
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

		zVec = _mm256_fmadd_pd(conc, chg, zVec);
		dZVec = _mm256_fmadd_pd(dConc, chg, dZVec);
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
#ifdef ECHMET_COMPILER_MSVC
template
void ChargeSummer<double, InstructionSet::FMA3, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR,
	const double *const ECHMET_RESTRICT_PTR,
	double &, double &) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
double ChargeSummer<double, InstructionSet::FMA3, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	VD vIs;
	__m256d isVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chgSq = M256D(m_chargesSquared + idx);

		isVec = _mm256_fmadd_pd(conc, chgSq, isVec);
	}

	_mm256_store_pd(vIs, isVec);

	is = vIs[0] + vIs[1] + vIs[2] + vIs[3];

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::FMA3, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

template <>
double ChargeSummer<double, InstructionSet::FMA3, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept
{
	double is;
	VD vIs;
	__m256d isVec = M256D(ZERO);

	size_t idx{0};
	for (; idx < m_NBlock; idx += m_blockSize) {
		__m256d conc = M256D(icConcs + idx);
		__m256d chgSq = M256D(m_chargesSquared + idx);

		isVec = _mm256_fmadd_pd(conc, chgSq, isVec);
	}

	_mm256_store_pd(vIs, isVec);

	is = vIs[0] + vIs[1] + vIs[2] + vIs[3];

	for (; idx < m_N; idx++)
		is += icConcs[idx] * m_chargesSquared[idx];

	return 0.0005 * is;
}
#ifdef ECHMET_COMPILER_MSVC
template
double ChargeSummer<double, InstructionSet::FMA3, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR) noexcept;
#endif // ECHMET_COMPILER_MSVC

} // namespace CAES
} // namespace ECHMET
