#ifndef ECHMET_MATH_NEWTON_RAPHSON_HPP
#define ECHMET_MATH_NEWTON_RAPHSON_HPP

namespace ECHMET {
namespace CAES {

template <typename NRReal>
void NewtonRaphson<NRReal>::AInit()
{}

template <typename NRReal>
NewtonRaphson<NRReal>::~NewtonRaphson()
{}

template <typename NRReal>
SolverVector<NRReal> const & NewtonRaphson<NRReal>::ASolve ()
{

	// INIT

	SolverVector<NRReal> &x = *m_px;

	m_iteration = 0;

	m_f = SolverVector<NRReal>::Zero(x.rows());
	m_j = SolverMatrix<NRReal>::Zero(x.rows(), x.cols());
	m_dx = SolverMatrix<NRReal>::Zero(x.rows(), x.cols());

	AInit();

	// DOIT

	while (true)
	{
		// !2
		ACalculateF(m_f, x);
		ACalculateJ(m_j, x);

		ZCalculateMeasures(m_f, m_fMin, m_fMax);
		ZCalculateMeasures(m_dx, m_dxMin, m_dxMax);                                  // 1

		ZCheckStatus();

		if (m_status != Status::CONTINUE)
			return x;

		this->m_iteration++;

		/* FIXME: Eigen does non throw if there is no solution!!! */
		try {
			//m_dx = arma::solve(m_j, m_f, arma::solve_opts::equilibrate + arma::solve_opts::no_approx);
			m_dx = m_j.lu().solve(m_f);
		} catch (std::runtime_error &) {
			m_status = Status::NO_SOLUTION;
			return x;
		}

		x -= m_dx;
	};

}
//  1 : it also applies to the 0th call when dxMin and dxMax only result in 0
// !2 : implementers may rely on the order of F/J evaluation
//---------------------------------------------------------------------------

template <typename NRReal>
NewtonRaphson<NRReal>::NewtonRaphson() :
	maxIterations(defaultMaxIterations()),
	xPrecision(defaultDxPrecision()),
	fPrecision(defaultFPrecision())
{
	ZConstructor();
}

template <typename NRReal>
NewtonRaphson<NRReal>::NewtonRaphson(SolverVector<NRReal> const &x_) :
	maxIterations(defaultMaxIterations()),
	xPrecision(defaultDxPrecision()),
	fPrecision(defaultFPrecision())
{
	ZConstructor();
	ZXAssign(x_);
}

template <typename NRReal>
NewtonRaphson<NRReal>::NewtonRaphson(SolverVector<NRReal> *x_) :
	maxIterations(defaultMaxIterations()),
	xPrecision(defaultDxPrecision()),
	fPrecision(defaultFPrecision())
{
	ZConstructor();
	m_px = x_;
}

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

template <typename NRReal>
void NewtonRaphson<NRReal>::ZXAssign(SolverVector<NRReal> const &ax)
{
	if (m_xInternal.size() != ax.size())
		m_xInternal.resize(ax.rows(), ax.cols());

	m_xInternal = ax;

	m_px = &m_xInternal;
}

} // namespace Math
} // namespace ECHMET

#endif // ECHMET_MATH_NEWTON_RAPHSON_HPP
