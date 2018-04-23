#ifndef ECHMET_CAES_SOLVERCONTEXTIMPL_H
#define ECHMET_CAES_SOLVERCONTEXTIMPL_H

#include <echmetcaes.h>
#include "types.h"
#include "totalequilibrium.h"

namespace ECHMET {
namespace CAES {

template <typename CAESReal>
class SolverContextImpl : public SolverContext {
public:
	explicit SolverContextImpl(const LigandVec<CAESReal> *allLigands, const LigandIonicFormVec<CAESReal> *allLigandIFs,
				   const CNVec<CAESReal> *complexNuclei, const FormVec<CAESReal> *allForms,
				   const SolverMatrix<CAESReal> *preJacobian,
				   const size_t concentrationCount,
				   const size_t analyticalConcentrationCount,
				   std::vector<TotalEquilibriumBase *> &&totalEquilibria,
				   const SolverContext::Options options,
				   const int TECount) noexcept;
	virtual ~SolverContextImpl() noexcept override;
	virtual void ECHMET_CC destroy() const noexcept override;
	virtual Options ECHMET_CC options() const noexcept override;
	virtual RetCode ECHMET_CC setOptions(const Options options) noexcept override;

	const LigandVec<CAESReal> *allLigands;			/*!< Vector of all ligands */
	const LigandIonicFormVec<CAESReal> *allLigandIFs;	/*!< Vector of all ligand ionic forms */
	const CNVec<CAESReal> *complexNuclei;			/*!< Vector of all complexing component */
	const FormVec<CAESReal> *allForms;			/*!< Vector of all generated complex forms */
	const SolverMatrix<CAESReal> *preJacobian;		/*!< Pregenerated part of the Jacobian */
	const size_t concentrationCount;			/*!< Total number of concentrations to be calculates */
	const size_t analyticalConcentrationCount;		/*!< Number of input analytical concentrations */
	const int TECount;					/*!< Number of ionic concentrations that can be estimated from G-polynomial */

	std::vector<TotalEquilibriumBase *> totalEquilibria;
	SolverVector<CAESReal> estimatedIonicConcentrations;
	SolverVector<CAESReal> dEstimatedIonicConcentrationsdH;

private:
	Options m_options;					/*!< Solver options */
};

} // namespace CAES
} // namespace ECHMET

#include "solvercontextimpl.hpp"

#endif // ECHMET_CAES_SOLVERCONTEXTIMPL_H
