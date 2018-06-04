#include "newtonraphson.h"

namespace ECHMET {
namespace CAES {

template <>
void NewtonRaphson<double, InstructionSet::FMA3>::initializeLu(const int elements)
{
	m_luCalc = new Eigen::PartialPivLU<SolverMatrix<double>>(elements);
}

template <>
typename NewtonRaphson<double, InstructionSet::FMA3>::TX const & NewtonRaphson<double, InstructionSet::FMA3>::ASolve()
{
	m_iteration = 0;
	m_stuckCounter = 0;

	AInit();

	// DOIT

	while (true)
	{
		// !2
		ACalculateF(m_f, *m_px);
		ACalculateJ(m_j, *m_px);

		ZCalculateMeasures(m_f, m_fMin, m_fMax);
		ZCalculateMeasures(m_dx, m_dxMin, m_dxMax);                                  // 1

		ZCheckStatus();

		if (m_status != Status::CONTINUE)
			return *m_px;

		this->m_iteration++;

		m_luCalc->compute(m_j);
		m_dx = m_luCalc->solve(m_f);

		*m_px -= m_dx;
	};
}

} // namespace CAES
} // namespace ECHMET
