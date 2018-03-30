#ifndef ECHMET_CAES_SOLVERINTERNAL_H
#define ECHMET_CAES_SOLVERINTERNAL_H

#include "solvercontextimpl.h"
#include "types.h"
#include "internals/newtonraphson.h"
#include <echmetcaes.h>

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
class SolverInternal : public NewtonRaphson<CAESReal> {
public:
	explicit SolverInternal(const SolverContextImpl<CAESReal> *ctx);
	const SolverContext * context() const noexcept;
	bool isReady() const noexcept;
	SolverIterations iterations() const noexcept;
	void makeReady();
	RetCode solve(const SolverVector<CAESReal> *analyticalConcentrations, const SolverVector<CAESReal> &estimatedConcentrations,  const bool isCorrection, const size_t iterations) noexcept;
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

	void ACalculateF(SolverVector<CAESReal> &Fx, const SolverVector<CAESReal> &pCx) override;
	void ACalculateJ(SolverMatrix<CAESReal> &Jx, const SolverVector<CAESReal> &pCx) override;
	CAESReal calculateIonicStrength() const;
	void initializepACoeffs();
	CAESReal pACoeff(const int charge);
	void recalculatepACoeffs(const CAESReal &is);
	void validateMatrix(const SolverMatrix<CAESReal> &m);
	void validateVector(const SolverVector<CAESReal> &v);

	const SolverContextImpl<CAESReal> *m_ctx;		/*!< Back pointer to the context used to initialize the object */
	const FormVec<CAESReal> *m_allForms;			/*!< Vector of all complex forms and free complex nuclei in the system */
	const LigandIonicFormVec<CAESReal> *m_allLigandIFs;	/*!< Vector of all ligand ionic forms */
	const LigandVec<CAESReal> *m_allLigands;		/*!< Vector of all ligands */
	const CNVec<CAESReal> *m_complexNuclei;			/*!< Vector of all complex nuclei */
	const SolverMatrix<CAESReal> *m_preJacobian;		/*!< Precalculated Jacobian */
	SolverVector<CAESReal> m_pCx;				/*!< Working vector with pX concentrations */
	SolverVector<CAESReal> m_rCx;				/*!< Working vector with X10 concentrations */
	const SolverVector<CAESReal> *m_analyticalConcentrations;
	CAESReal m_finalIonicStrength;				/*!< Ionic strength of the resolved system */

	uint32_t m_outerIterations;
	uint32_t m_totalIterations;

	bool m_correctForIS;					/*!< Ionic strength correction On/Off */
	std::vector<CAESReal> m_pACoeffs;			/*!< Activity coefficients for all charges */
	std::vector<int> m_chargesSquared;			/*!< Vector of squared charges */
	int m_chargeMax;					/*!< Absolute value of maximum charge present in the system */
};


} // namespace CAES
} // namespace ECHMET

#include "solverinternal.hpp"

#endif // ECHMET_CAES_SOLVERINTERNAL_H


