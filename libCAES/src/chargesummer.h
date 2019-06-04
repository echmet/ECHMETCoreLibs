#ifndef ECHMET_CAES_CHARGESUMMER_H
#define ECHMET_CAES_CHARGESUMMER_H

#include "solvercontextimpl.h"
#include "totalequilibrium.h"
#include "mappedmatrix.h"

namespace ECHMET {
namespace CAES {

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
class ChargeSummer {
public:
	using MappedVector = typename MMTypes<CAESReal, ISet>::Vector;

	ChargeSummer(const size_t N, const std::vector<TotalEquilibriumBase *> &totalEquilibria) :
		m_N{N},
		m_blockSize{VDType<ISet>::ALIGNMENT_BYTES / sizeof(CAESReal)},
		m_NBlock{N - (N % m_blockSize)}
	{
		m_charges = AlignedAllocator<double, VDType<ISet>::ALIGNMENT_BYTES>::alloc(m_N);
		m_chargesSquared = AlignedAllocator<double, VDType<ISet>::ALIGNMENT_BYTES>::alloc(m_N);

		size_t ctr{2};
		for (const TotalEquilibriumBase *teb : totalEquilibria) {
			const auto *te = static_cast<const TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
			for (int charge = te->numLow; charge <= te->numHigh; charge++) {
				m_charges[ctr] = charge;
				m_chargesSquared[ctr] = charge * charge;

				ctr++;
			}
		}
		m_charges[0] = 1.0;	/* H+ */
		m_charges[1] = -1.0;	/* OH- */

		m_chargesSquared[0] = 1.0;	/* H+ */
		m_chargesSquared[1] = 1.0;	/* OH- */
	}


	~ChargeSummer() noexcept
	{
		AlignedAllocator<double, VDType<ISet>::ALIGNMENT_BYTES>::free(m_charges);
		AlignedAllocator<double, VDType<ISet>::ALIGNMENT_BYTES>::free(m_chargesSquared);
	}

	CAESReal calc(const CAESReal *const ECHMET_RESTRICT_PTR icConcs) noexcept
	{
		CAESReal z{0};

		for (size_t idx = 0; idx < m_N; idx++)
			z += m_charges[idx] * icConcs[idx];

		return z;
	}

	CAESReal calculateIonicStrength(const CAESReal *const ECHMET_RESTRICT_PTR icConcs) noexcept
	{
		CAESReal ionicStrength{0};

		for (size_t idx = 0; idx < m_N; idx++)
			ionicStrength += m_chargesSquared[idx] * icConcs[idx];

		return 0.0005 * ionicStrength; /* Scale to mol/dm3 */
	}

	void calcWithdZ(const CAESReal *const ECHMET_RESTRICT_PTR icConcs,
			const CAESReal *const ECHMET_RESTRICT_PTR dIcConcsdH,
			CAESReal &z, CAESReal &dZ) noexcept
	{
		z = CAESReal{0};
		dZ = CAESReal{0};

		for (size_t idx = 0; idx < m_N; idx++) {
			z += m_charges[idx] * icConcs[idx];
			dZ += m_charges[idx] * dIcConcsdH[idx];
		}
	}

private:
	const size_t m_N;
	const size_t m_blockSize;
	const size_t m_NBlock;

	double *m_charges;
	double *m_chargesSquared;
};

template <>
double ChargeSummer<double, InstructionSet::SSE2, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::SSE2, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
void ChargeSummer<double, InstructionSet::SSE2, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								   const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								   double &z, double &dZ) noexcept;
template <>
void ChargeSummer<double, InstructionSet::SSE2, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept;
template <>
double ChargeSummer<double, InstructionSet::SSE2, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::SSE2, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;

template <>
double ChargeSummer<double, InstructionSet::AVX, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::AVX, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
void ChargeSummer<double, InstructionSet::AVX, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept;
template <>
void ChargeSummer<double, InstructionSet::AVX, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								 const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								 double &z, double &dZ) noexcept;

template <>
double ChargeSummer<double, InstructionSet::FMA3, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::FMA3, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
void ChargeSummer<double, InstructionSet::FMA3, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								   const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								   double &z, double &dZ) noexcept;
template <>
void ChargeSummer<double, InstructionSet::FMA3, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept;

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CHARGESUMMER_H
