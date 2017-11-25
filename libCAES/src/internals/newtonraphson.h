#ifndef ECHMET_MATH_NEWTON_RAPHSON_H
#define ECHMET_MATH_NEWTON_RAPHSON_H

#include <cstddef>
#include <internal/echmetmath_internal.h>
#include "../types.h"

namespace ECHMET {
namespace CAES {

/* Prevent multiple definitions during linking */
namespace {

template <typename NRReal>
NRReal getDefaultDxPrecision()
{
	return 1.0e-12;
}

template <typename NRReal>
NRReal getDefaultFPrecision()
{
	return 1.0e-12;
}

template <>
mpfr::mpreal getDefaultDxPrecision()
{
	return "1.0e-50";
}

template <>
mpfr::mpreal getDefaultFPrecision()
{
	return "1.0e-50";
}

}

template <typename NRReal>
class NewtonRaphson {
public:
	enum class Status {
		NONE,
		CONTINUE,
		SUCCEEDED,
		NO_CONVERGENCE,
		STUCK,
		NO_SOLUTION,
		ABORTED
	};

	size_t maxIterations;
	NRReal const xPrecision;
	NRReal const fPrecision;

	NewtonRaphson();
	NewtonRaphson(SolverVector<NRReal> const & matrix);
	NewtonRaphson(SolverVector<NRReal> *matrix);
	NewtonRaphson(NewtonRaphson const &) = delete;
	virtual ~NewtonRaphson();

	inline void               Report()        const;

	//SolverMatrix const &  GetX()    const;

	int                GetIterations() const;
	Status             GetStatus()     const;
	//SolverMatrix const &  GetF()          const;
	SolverMatrix<NRReal> const &  GetJ()          const;
	SolverMatrix<NRReal> const &  GetdX()         const;

	NRReal             GetMaxdx()      const;
	NRReal             GetMindx()      const;
	NRReal             GetMaxF()       const;
	NRReal             GetMinF()       const;

	void operator=(NewtonRaphson const &) = delete;

private:
	Status m_status;
	NRReal m_dxMin;
	NRReal m_dxMax;
	NRReal m_fMin;
	NRReal m_fMax;

	SolverVector<NRReal> *m_px;

	SolverVector<NRReal> m_xInternal;
	SolverVector<NRReal> m_f;
	SolverMatrix<NRReal> m_j;
	SolverMatrix<NRReal> m_dx;

	void ZConstructor();
	void ZXAssign(SolverVector<NRReal> const &ax);
	void ZCalculateMeasures(SolverMatrix<NRReal> const &m, NRReal &min, NRReal &max);
	void ZCheckStatus();

	static size_t defaultMaxIterations() noexcept
	{
		return 1000;
	}

	static NRReal defaultDxPrecision() noexcept
	{
		return getDefaultDxPrecision<NRReal>();
	}

	static NRReal defaultFPrecision() noexcept
	{
		return getDefaultFPrecision<NRReal>();
	}

protected:
	virtual void AInit();
	virtual void ACalculateF(SolverVector<NRReal> &, SolverVector<NRReal> const &) = 0;
	virtual void ACalculateJ(SolverMatrix<NRReal> &, SolverVector<NRReal> const &) = 0;
	SolverVector<NRReal> const & ASolve();
	inline SolverVector<NRReal> const & ASolve(SolverVector<NRReal> const &);
	inline SolverVector<NRReal> const & ASolve(SolverVector<NRReal> *);
	size_t m_iteration;
};

template <typename NRReal>
void NewtonRaphson<NRReal>::ZConstructor()
{
	m_iteration = 0;
	m_status    = Status::NONE;
	m_dxMin     = 0.0;
	m_dxMax     = 0.0;
	m_fMin      = 0.0;
	m_fMax      = 0.0;

	m_px          = nullptr; // 1
}
//---------------------------------------------------------------------------
// 1) We may set $px to $x_internal. But setting it to NULL leads to an exception
//    if the solver is not initialized at all.

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
template <typename NRReal>
void NewtonRaphson<NRReal>::ZCalculateMeasures(SolverMatrix<NRReal> const &m, NRReal &min, NRReal &max)
{
	min = VMath::abs(m(0, 0));
	max = min;

	for (int row = 1; row < m.rows(); row++) {
		const SolverVector<NRReal> &v = m.row(row);

		const NRReal val = VMath::abs(v(0));

		if (val > max)
			max = val;
		else if (val < min)
			min = val;
	}
}

template <typename NRReal>
void NewtonRaphson<NRReal>::ZCheckStatus()
{
	// order of ifs matters

	if (m_fMax <= fPrecision && m_dxMax <= xPrecision) m_status = Status::SUCCEEDED;
	else if (m_iteration >= maxIterations)             m_status = Status::NO_CONVERGENCE;
	else if (m_iteration == 0)                         m_status = Status::CONTINUE;
	else if (m_dxMax <= xPrecision)                    m_status = Status::STUCK;
	else                                               m_status = Status::CONTINUE;
}

template <typename NRReal>
SolverVector<NRReal> const & NewtonRaphson<NRReal>::ASolve(SolverVector<NRReal> const &x_)
{
	ZXAssign(x_);

	return ASolve();
}

template <typename NRReal>
SolverVector<NRReal> const & NewtonRaphson<NRReal>::ASolve(SolverVector<NRReal> *x_)
{
	m_px = x_;

	return ASolve();
}


} // namespace MATH
} // namespace ECHMET

#include "newtonraphson.hpp"

#endif // ECHMET_MATH_NEWTON_RAPHSON_H
