#ifndef ECHMET_CAES_CHARGESUMMER_H
#define ECHMET_CAES_CHARGESUMMER_H

#include "solvercontextimpl.h"
#include "totalequilibrium.h"
#ifdef ECHMET_USE_X86_EXTENSIONS
#include "vecmath/vecmath.h"
#else
#include "genericmath.h"
#endif // ECHMET_USE_X86_EXTENSIONS

#include <cassert>

#include <iostream>

namespace ECHMET {
namespace CAES {

class BlockSize {
public:
	template <typename CAESReal, InstructionSet ISet, std::enable_if_t<std::is_fundamental<CAESReal>::value, int> = 0>
	static
	size_t blockSize() noexcept
	{
		return VDType<ISet>::ALIGNMENT_BYTES / sizeof(CAESReal);
	}

	template <typename CAESReal, InstructionSet ISet, std::enable_if_t<!std::is_fundamental<CAESReal>::value, int> = 0>
	static
	size_t blockSize() noexcept
	{
		return 0;
	}

	template <typename CAESReal, std::enable_if_t<std::is_fundamental<CAESReal>::value, int> = 0>
	static
	size_t nBlock(const size_t N, const size_t blockSize) noexcept
	{
		assert(blockSize > 0);

		return N - (N % blockSize);
	}

	template <typename CAESReal, std::enable_if_t<!std::is_fundamental<CAESReal>::value, int> = 0>
	static
	size_t nBlock(const size_t, const size_t) noexcept
	{
		return 0;
	}
};

template <typename CAESReal, InstructionSet ISet, bool ThreadSafe>
class ChargeSummer {
public:
	ChargeSummer(const size_t N, const std::vector<TotalEquilibrium<CAESReal, ISet, ThreadSafe>> &totalEquilibria) :
		m_N{N},
		m_blockSize{BlockSize::blockSize<CAESReal, ISet>()},
		m_NBlock{BlockSize::nBlock<CAESReal>(m_N, m_blockSize)}
	{
		m_charges = AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(m_N);
		m_chargesSquared = AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(m_N);

		size_t ctr{2};
		for (const auto &te : totalEquilibria) {
			for (int charge = te.numLow; charge <= te.numHigh; charge++) {
				m_charges[ctr] = charge;
				m_chargesSquared[ctr] = charge * charge;

				ctr++;
			}
		}
		m_charges[0] = CAESReal{1.0};	/* H+ */
		m_charges[1] = CAESReal{-1.0};	/* OH- */

		m_chargesSquared[0] = CAESReal{1.0};	/* H+ */
		m_chargesSquared[1] = CAESReal{1.0};	/* OH- */
	}

	ChargeSummer(const ChargeSummer &) = delete;
	ChargeSummer(const ChargeSummer &&) = delete;

	~ChargeSummer() noexcept
	{
		AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_charges);
		AlignedAllocator<CAESReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_chargesSquared);
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

	CAESReal *ECHMET_RESTRICT_PTR m_charges;
	CAESReal *ECHMET_RESTRICT_PTR m_chargesSquared;
};

#ifdef ECHMET_USE_X86_EXTENSIONS
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
double ChargeSummer<double, InstructionSet::SSE2, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;

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
double ChargeSummer<double, InstructionSet::AVX, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::AVX, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;

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
template <>
double ChargeSummer<double, InstructionSet::FMA3, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::FMA3, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;

#ifndef ECHMET_DISABLE_AVX512

template <>
double ChargeSummer<double, InstructionSet::AVX512, false>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::AVX512, true>::calc(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
void ChargeSummer<double, InstructionSet::AVX512, false>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								   const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								   double &z, double &dZ) noexcept;
template <>
void ChargeSummer<double, InstructionSet::AVX512, true>::calcWithdZ(const double *const ECHMET_RESTRICT_PTR icConcs,
								  const double *const ECHMET_RESTRICT_PTR dIcConcsdH,
								  double &z, double &dZ) noexcept;
template <>
double ChargeSummer<double, InstructionSet::AVX512, false>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;
template <>
double ChargeSummer<double, InstructionSet::AVX512, true>::calculateIonicStrength(const double *const ECHMET_RESTRICT_PTR icConcs) noexcept;

#endif // ECHMET_DISABLE_AVX512

#endif // ECHMET_USE_X86_EXTENSIONS

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CHARGESUMMER_H
