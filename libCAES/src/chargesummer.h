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
		m_chargesArray = AlignedAllocator<int, VDType<ISet>::ALIGNMENT_BYTES>::alloc(m_N);

		size_t ctr{2};
		for (const TotalEquilibriumBase *teb : totalEquilibria) {
			const auto *te = static_cast<const TotalEquilibrium<CAESReal, ThreadSafe> *>(teb);
			for (int charge = te->numLow; charge <= te->numHigh; charge++)
				m_chargesArray[ctr++] = charge;

		}
		m_chargesArray[0] = 1;	/* H+ */
		m_chargesArray[1] = -1;	/* OH- */
	}


	~ChargeSummer() noexcept
	{
		AlignedAllocator<int, VDType<ISet>::ALIGNMENT_BYTES>::free(m_chargesArray);
	}

	CAESReal calc(const MappedVector &icConcs) noexcept
	{
		CAESReal z{0};

		for (size_t idx = 0; idx < m_N; idx++)
			z += m_chargesArray[idx] * icConcs(idx);

		return z;
	}

	void calcWithdZ(const MappedVector &icConcs, const MappedVector &dIcConcsdH,
			CAESReal &z, CAESReal &dZ) noexcept
	{
		z = CAESReal{0};
		dZ = CAESReal{0};

		for (size_t idx = 0; idx < m_N; idx++) {
			z += m_chargesArray[idx] * icConcs(idx);
			dZ += m_chargesArray[idx] * dIcConcsdH(idx);
		}
	}

private:
	const size_t m_N;
	const size_t m_blockSize;
	const size_t m_NBlock;

	int *m_chargesArray;
};

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_CHARGESUMMER_H
