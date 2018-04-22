#ifndef ECHMET_MATH_NEWTON_RAPHSON_HPP
#define ECHMET_MATH_NEWTON_RAPHSON_HPP

namespace ECHMET {
namespace CAES {

template <typename NRReal>
NewtonRaphson<NRReal>::NewtonRaphson(const int elements) :
	maxIterations(defaultMaxIterations()),
	xPrecision(defaultDxPrecision()),
	fPrecision(defaultFPrecision()),
	m_lu(elements)
{
	ZConstructor();

	m_f = SolverVector<NRReal>::Zero(elements);
	m_j = SolverMatrix<NRReal>::Zero(elements, elements);
	m_dx = SolverVector<NRReal>::Zero(elements);

	ECHMET_DEBUG_CODE(fprintf(stderr, "Fx rows %ld, Jx dims %ld,%ld, Dx rows %ld\n", m_f.rows(), m_j.rows(), m_j.cols(), m_dx.rows()));
}

template <typename NRReal>
NewtonRaphson<NRReal>::~NewtonRaphson()
{}

template <typename NRReal>
void NewtonRaphson<NRReal>::AInit()
{}

template <typename NRReal>
typename NewtonRaphson<NRReal>::TX const & NewtonRaphson<NRReal>::ASolve()
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

		//m_dx = m_j.partialPivLu().solve(m_f);
		m_lu.compute(m_j);
		m_dx = m_lu.solve(m_f);

		*m_px -= m_dx;
	};

}
//  1 : it also applies to the 0th call when dxMin and dxMax only result in 0
// !2 : implementers may rely on the order of F/J evaluation
//---------------------------------------------------------------------------

template <typename NRReal>
int NewtonRaphson<NRReal>::GetIterations() const
{ return m_iteration; }

template <typename NRReal>
typename NewtonRaphson<NRReal>::Status NewtonRaphson<NRReal>::GetStatus() const
{ return m_status; }

/*SolverMatrix const &  NewtonRaphson::GetF() const
{ return m_f; }*/

template <typename NRReal>
SolverMatrix<NRReal> const & NewtonRaphson<NRReal>::GetJ() const
{ return m_j; }

template <typename NRReal>
SolverMatrix<NRReal> const & NewtonRaphson<NRReal>::GetdX() const
{ return m_dx; }

template <typename NRReal>
NRReal NewtonRaphson<NRReal>::GetMindx() const
{ return m_dxMin; }

template <typename NRReal>
NRReal NewtonRaphson<NRReal>::GetMaxdx() const
{ return m_dxMax; }

template <typename NRReal>
NRReal NewtonRaphson<NRReal>::GetMinF() const
{ return m_fMin; }

template <typename NRReal>
NRReal NewtonRaphson<NRReal>::GetMaxF() const
{ return m_fMax; }

} // namespace Math
} // namespace ECHMET

#endif // ECHMET_MATH_NEWTON_RAPHSON_HPP
