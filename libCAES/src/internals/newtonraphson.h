#ifndef ECHMET_MATH_NEWTON_RAPHSON_H
#define ECHMET_MATH_NEWTON_RAPHSON_H

#include <cstddef>
#include <internal/echmetmath_internal.h>
#include "../types.h"
#include "../vecmath/vecmath.h"

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

template <typename NRReal, InstructionSet ISet>
class NRTypes {
public:
	typedef Eigen::Map<SolverMatrix<NRReal>, Eigen::AlignmentType::Aligned16> Matrix;
	typedef Eigen::Map<SolverVector<NRReal>, Eigen::AlignmentType::Aligned16> Vector;
};

template <>
class NRTypes<double, InstructionSet::AVX> {
public:
	typedef Eigen::Map<SolverMatrix<double>, Eigen::AlignmentType::Aligned32> Matrix;
	typedef Eigen::Map<SolverVector<double>, Eigen::AlignmentType::Aligned32> Vector;
};

template <>
class NRTypes<double, InstructionSet::FMA3> {
public:
	typedef Eigen::Map<SolverMatrix<double>, Eigen::AlignmentType::Aligned32> Matrix;
	typedef Eigen::Map<SolverVector<double>, Eigen::AlignmentType::Aligned32> Vector;
};

template <typename NRReal, InstructionSet ISet>
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

	typedef typename NRTypes<NRReal, ISet>::Matrix TM;
	typedef typename NRTypes<NRReal, ISet>::Vector TX;

	size_t maxIterations;
	NRReal const xPrecision;
	NRReal const fPrecision;

	explicit NewtonRaphson(const int elements);
	NewtonRaphson(NewtonRaphson const &) = delete;
	virtual ~NewtonRaphson();

	inline void        Report()        const;

	int                GetIterations() const;
	Status             GetStatus()     const;
	TM const &         GetJ()          const;
	TX const &         GetdX()         const;

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

	TX *m_px;

	NRReal *m_f_raw;
	NRReal *m_j_raw;
	NRReal *m_dx_raw;

	TX m_f;
	TM m_j;
	TX m_dx;

	Eigen::PartialPivLU<SolverMatrix<NRReal>> *m_luCalc;

	int m_stuckCounter;

	void initializeLu(const int elements);

	void ZConstructor();
	void ZCalculateMeasures(TX const &v, NRReal &min, NRReal &max);
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
	virtual void ACalculateF(TX &, TX const &) = 0;
	virtual void ACalculateJ(TM &, TX const &) = 0;
	TX const & ASolve();
	TX const & ASolve(TX *);
	size_t m_iteration;
};

template <typename NRReal, InstructionSet ISet>
void NewtonRaphson<NRReal, ISet>::ZConstructor()
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
template <typename NRReal, InstructionSet ISet>
void NewtonRaphson<NRReal, ISet>::ZCalculateMeasures(TX const &v, NRReal &min, NRReal &max)
{
	min = VMath::abs(v(0));
	max = min;

	for (int row = 1; row < v.rows(); row++) {
		const NRReal val = VMath::abs(v(row));

		if (val > max)
			max = val;
		else if (val < min)
			min = val;
	}
}

template <typename NRReal, InstructionSet ISet>
void NewtonRaphson<NRReal, ISet>::ZCheckStatus()
{
	// order of ifs matters

	if (m_fMax <= fPrecision /*&& m_dxMax <= xPrecision*/)
		m_status = Status::SUCCEEDED;
	else if (m_iteration >= maxIterations)
		m_status = Status::NO_CONVERGENCE;
	else if (m_iteration == 0)
		m_status = Status::CONTINUE;
	else if (m_dxMax <= xPrecision && m_stuckCounter++ > 20)
		m_status = Status::STUCK;
	else {
		m_status = Status::CONTINUE;
		m_stuckCounter = 0;
	}
}

template <typename NRReal, InstructionSet ISet>
typename NewtonRaphson<NRReal, ISet>::TX const & NewtonRaphson<NRReal, ISet>::ASolve(TX *x_)
{
	m_px = x_;

	return ASolve();
}

} // namespace MATH
} // namespace ECHMET

#include "newtonraphson.hpp"

#endif // ECHMET_MATH_NEWTON_RAPHSON_H
