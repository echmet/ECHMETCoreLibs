#ifndef ECHMET_MATH_NEWTON_RAPHSON_HPP
#define ECHMET_MATH_NEWTON_RAPHSON_HPP

namespace ECHMET {
namespace CAES {

#ifdef ECHMET_USE_X86_EXTENSIONS
template <>
void NewtonRaphson<double, InstructionSet::AVX>::initializeLu(const int elements);

template <>
void NewtonRaphson<double, InstructionSet::FMA3>::initializeLu(const int elements);

#ifndef ECHMET_DISABLE_AVX512

template <>
void NewtonRaphson<double, InstructionSet::AVX512>::initializeLu(const int elements);

#endif // ECHMET_DISABLE_AVX512

#endif // ECHMET_USE_X86_EXTENSIONS

template <typename NRReal, InstructionSet ISet>
void NewtonRaphson<NRReal, ISet>::initializeLu(const int elements)
{
	m_luCalc = new Eigen::PartialPivLU<SolverMatrix<NRReal>>(elements);
}

template <typename NRReal, InstructionSet ISet>
NewtonRaphson<NRReal, ISet>::NewtonRaphson(const int elements) :
	maxIterations(defaultMaxIterations()),
	xPrecision(defaultDxPrecision()),
	fPrecision(defaultFPrecision()),
	m_f_raw(AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(elements)),
	m_j_raw(AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(elements * elements)),
	m_dx_raw(AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::alloc(elements)),
	m_f(m_f_raw, elements),
	m_j(m_j_raw, elements, elements),
	m_dx(m_dx_raw, elements)
{
	if (m_f_raw == nullptr || m_j_raw == nullptr || m_dx_raw == nullptr) {
		AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_f_raw);
		AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_j_raw);
		AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_dx_raw);

		throw std::bad_alloc();
	}

	ZConstructor();

	try {
		initializeLu(elements);
	} catch (std::bad_alloc &) {
		AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_f_raw);
		AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_j_raw);
		AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_dx_raw);

		throw;
	}

	ECHMET_DEBUG_CODE(fprintf(stderr, "Fx rows %ld, Jx dims %ld,%ld, Dx rows %ld\n", m_f.rows(), m_j.rows(), m_j.cols(), m_dx.rows()));
}

template <typename NRReal, InstructionSet ISet>
NewtonRaphson<NRReal, ISet>::~NewtonRaphson()
{
	AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_f_raw);
	AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_j_raw);
	AlignedAllocator<NRReal, VDType<ISet>::ALIGNMENT_BYTES>::free(m_dx_raw);

	delete m_luCalc;
}

template <typename NRReal, InstructionSet ISet>
void NewtonRaphson<NRReal, ISet>::AInit()
{}

#ifdef ECHMET_USE_X86_EXTENSIONS
template <>
typename NewtonRaphson<double, InstructionSet::AVX>::TX const & NewtonRaphson<double, InstructionSet::AVX>::ASolve();

template <>
typename NewtonRaphson<double, InstructionSet::FMA3>::TX const & NewtonRaphson<double, InstructionSet::FMA3>::ASolve();

template <>
typename NewtonRaphson<double, InstructionSet::AVX512>::TX const & NewtonRaphson<double, InstructionSet::AVX512>::ASolve();
#endif // ECHMET_USE_X86_EXTENSIONS

template <typename NRReal, InstructionSet ISet>
typename NewtonRaphson<NRReal, ISet>::TX const & NewtonRaphson<NRReal, ISet>::ASolve()
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

		m_luCalc->compute(m_j);
		m_dx = m_luCalc->solve(m_f);

		ZCalculateMeasures(m_f, m_fMin, m_fMax);
		ZCalculateMeasures(m_dx, m_dxMin, m_dxMax);

		ZCheckStatus();

		if (m_status != Status::CONTINUE)
			return *m_px;

		this->m_iteration++;

		*m_px -= m_dx;
	};
}

// !2 : implementers may rely on the order of F/J evaluation
//---------------------------------------------------------------------------

template <typename NRReal, InstructionSet ISet>
int NewtonRaphson<NRReal, ISet>::GetIterations() const
{ return m_iteration; }

template <typename NRReal, InstructionSet ISet>
typename NewtonRaphson<NRReal, ISet>::Status NewtonRaphson<NRReal, ISet>::GetStatus() const
{ return m_status; }

/*SolverMatrix const &  NewtonRaphson::GetF() const
{ return m_f; }*/

template <typename NRReal, InstructionSet ISet>
typename NewtonRaphson<NRReal, ISet>::TM const & NewtonRaphson<NRReal, ISet>::GetJ() const
{ return m_j; }

template <typename NRReal, InstructionSet ISet>
typename NewtonRaphson<NRReal, ISet>::TX const & NewtonRaphson<NRReal, ISet>::GetdX() const
{ return m_dx; }

template <typename NRReal, InstructionSet ISet>
NRReal NewtonRaphson<NRReal, ISet>::GetMindx() const
{ return m_dxMin; }

template <typename NRReal, InstructionSet ISet>
NRReal NewtonRaphson<NRReal, ISet>::GetMaxdx() const
{ return m_dxMax; }

template <typename NRReal, InstructionSet ISet>
NRReal NewtonRaphson<NRReal, ISet>::GetMinF() const
{ return m_fMin; }

template <typename NRReal, InstructionSet ISet>
NRReal NewtonRaphson<NRReal, ISet>::GetMaxF() const
{ return m_fMax; }

} // namespace Math
} // namespace ECHMET

#endif // ECHMET_MATH_NEWTON_RAPHSON_HPP
