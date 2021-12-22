#ifndef ECHMET_CAES_SOLVERINTERNAL_H
#define ECHMET_CAES_SOLVERINTERNAL_H

#include "solvercontextimpl.h"
#include "internals/newtonraphson.h"
#ifdef ECHMET_USE_X86_EXTENSIONS
#include "vecmath/vecmath.h"
#else
#include "genericmath.h"
#endif // ECHMET_USE_X86_EXTENSIONS
#include <echmetcaes.h>

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
class NumericPrecisionSetter {};

template <>
class NumericPrecisionSetter<mpfr::mpreal> {
public:
	NumericPrecisionSetter()
	{
		m_currentPrecision = mpfr::mpreal::get_default_prec();

		if (m_currentPrecision != BITS_OF_PREC)
			mpfr::mpreal::set_default_prec(BITS_OF_PREC);
	}

	~NumericPrecisionSetter()
	{
		mpfr::mpreal::set_default_prec(m_currentPrecision);
	}

	static const long BITS_OF_PREC = 665; /* Corresponds to 200 digits of precision */

private:
	int m_currentPrecision;
};

template <typename CAESReal, InstructionSet ISet>
class SolverInternal : public NewtonRaphson<CAESReal, ISet> {
public:
	explicit SolverInternal(const SolverContextImpl<CAESReal> *ctx);
	~SolverInternal();
	const SolverContext * context() const noexcept;
	SolverIterations iterations() const noexcept;
	RetCode solve(const SolverVector<CAESReal> *analyticalConcentrations, const CAESReal *estimatedConcentrations,
		      const bool isCorrection, const size_t iterations, const CAESReal &ionicStrength) noexcept;
	SolverVector<CAESReal> rawConcentrations() const;
	CAESReal rawIonicStrength() const;
	void resultsToOutput(SysComp::CalculatedProperties &calcProps) noexcept;

private:
	class NumericErrorException : public std::exception {
	public:
		enum class ExType {
			ET_INF,
			ET_NAN
		};

		NumericErrorException(const ExType t);
		const char * what() const noexcept;

	private:
		const ExType m_t;

		static const char *infError;			/*<! Infinity error string */
		static const char *nanError;			/*<! Not-a-number error string */
	};

	template <InstructionSet MISet>
	class VectorizedDelogifier {
	public:
		explicit VectorizedDelogifier(const VecMath<MISet> &vecMath, const size_t N) :
			m_vecMath(vecMath),
			N(N),
			blockSize(calcBlockSize()),
			NBlock(N - (N % blockSize))
		{
		}

		void operator()(CAESReal *ECHMET_RESTRICT_PTR dst, const CAESReal *ECHMET_RESTRICT_PTR src)
		{
			ECHMET_DEBUG_CODE(fprintf(stderr, "Generic delogifier\n"));

			for (size_t idx = 0; idx < N; idx++)
				dst[idx] = X10(src[idx]);
		}

	private:
		size_t calcBlockSize()
		{
			const size_t sz = VDType<MISet>::ALIGNMENT_BYTES / sizeof(CAESReal);

			return sz > 0 ? sz : 1;
		}

		const VecMath<MISet> &m_vecMath;
		const size_t N;
		const size_t blockSize;
		const size_t NBlock;
	};

	template <InstructionSet MISet>
	class VectorizedLogifier {
	public:
		explicit VectorizedLogifier(const VecMath<MISet> &vecMath, const size_t N) :
			m_vecMath(vecMath),
			N(N),
			blockSize(VDType<MISet>::ALIGNMENT_BYTES / sizeof(CAESReal)),
			NBlock(N - (N % blockSize))
		{}

		void operator()(CAESReal *ECHMET_RESTRICT_PTR dst, const CAESReal *ECHMET_RESTRICT_PTR src)
		{
			ECHMET_DEBUG_CODE(fprintf(stderr, "Generic logifier\n"));

			for (size_t idx = 0; idx < N; idx++)
				dst[idx] = pX(src[idx]);
		}

	private:
		const VecMath<MISet> &m_vecMath;
		const size_t N;
		const size_t blockSize;
		const size_t NBlock;
	};

	void ACalculateF(typename NewtonRaphson<CAESReal, ISet>::TX &Fx, const typename NewtonRaphson<CAESReal, ISet>::TX &pCx) override;
	void ACalculateJ(typename NewtonRaphson<CAESReal, ISet>::TM &Jx, const typename NewtonRaphson<CAESReal, ISet>::TX &pCx) override;
	CAESReal calculateIonicStrength() const;
	void initializepACoeffs();
	CAESReal pACoeff(const int charge);
	void recalculatepACoeffs(const CAESReal &is);
	void validateMatrix(const SolverMatrix<CAESReal> &m);
	void validateVector(const typename NewtonRaphson<CAESReal, ISet>::TX &v);
	void validateVector(const SolverVector<CAESReal> &v);

	const SolverContextImpl<CAESReal> *m_ctx;			/*!< Back pointer to the context used to initialize the object */
	const FormVec<CAESReal> *m_allForms;				/*!< Vector of all complex forms and free complex nuclei in the system */
	const LigandIonicFormVec<CAESReal> *m_allLigandIFs;		/*!< Vector of all ligand ionic forms */
	const LigandVec<CAESReal> *m_allLigands;			/*!< Vector of all ligands */
	const CNVec<CAESReal> *m_complexNuclei;				/*!< Vector of all complex nuclei */
	const SolverMatrix<CAESReal> *m_preJacobian;			/*!< Precalculated Jacobian */
	const SolverVector<CAESReal> *m_analyticalConcentrations;	/*!< Pointer to vector of input analytical concentrations */
	CAESReal m_finalIonicStrength;					/*!< Ionic strength of the resolved system */

	CAESReal *m_pCx_raw;						/*!< Raw memory with pCx vector data */
	CAESReal *m_rCx_raw;						/*!< Raw memory with rCx vector data */
	typename NewtonRaphson<CAESReal, ISet>::TX m_pCx;		/*!< Working vector with pX concentrations */
	typename NewtonRaphson<CAESReal, ISet>::TX m_rCx;		/*!< Working vector with X10 concentrations */

	VecMath<ISet> *m_vecMath;
	VectorizedDelogifier<ISet> m_vecDelog;
	//VectorizedLogifier<ISet> m_vecLog;	<- Currently unused */

	uint32_t m_outerIterations;
	uint32_t m_totalIterations;

	bool m_correctForIS;					/*!< Ionic strength correction On/Off */
	std::vector<CAESReal> m_pACoeffs;			/*!< Activity coefficients for all charges */
};

} // namespace CAES
} // namespace ECHMET

#include "solverinternal.hpp"

#endif // ECHMET_CAES_SOLVERINTERNAL_H


