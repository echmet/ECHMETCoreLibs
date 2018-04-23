#ifndef ECHMET_CAES_SOLVERCONTEXTIMPL_HPP
#define ECHMET_CAES_SOLVERCONTEXTIMPL_HPP

#include "funcs.h"

namespace ECHMET {
namespace CAES {

/*!
 * SolverContext c-tor.
 *
 * @param[in] allLigands Pointer to vector of all ligands in the system.
 * @param[in] allLigandIFs Pointer to vector of all ligand ionic forms in the system.
 * @param[in] complexNuclei Pointer to vector of all complex nuclei in the system.
 * @param[in] allForms Pointer to all complex forms in the system.
 * @param[in] preJacobian Pointer to precomputed Jacobian.
 */
template <typename CAESReal>
SolverContextImpl<CAESReal>::SolverContextImpl(const LigandVec<CAESReal> *allLigands, const LigandIonicFormVec<CAESReal> *allLigandIFs,
					       const CNVec<CAESReal> *complexNuclei, const FormVec<CAESReal> *allForms,
					       const SolverMatrix<CAESReal> *preJacobian,
					       const size_t concentrationCount,
					       const size_t analyticalConcentrationCount,
					       std::vector<TotalEquilibriumBase *> &&totalEquilibria,
					       const SolverContext::Options options,
					       const int TECount) noexcept :
	allLigands(allLigands),
	allLigandIFs(allLigandIFs),
	complexNuclei(complexNuclei),
	allForms(allForms),
	preJacobian(preJacobian),
	concentrationCount(concentrationCount),
	analyticalConcentrationCount(analyticalConcentrationCount),
	TECount(TECount),
	totalEquilibria(totalEquilibria),
	m_options(options)
{
	estimatedIonicConcentrations.resize(TECount);
	dEstimatedIonicConcentrationsdH.resize(TECount);
}

/*!
 * SolverContext d-tor.
 */
template <typename CAESReal>
SolverContextImpl<CAESReal>::~SolverContextImpl() noexcept
{
	releasePointerContainer(allLigands);
	releasePointerContainer(allLigandIFs);
	releasePointerContainer(complexNuclei);
	releasePointerContainer(allForms);

	delete preJacobian;
}

template <typename CAESReal>
void ECHMET_CC SolverContextImpl<CAESReal>::destroy() const noexcept
{
	delete this;
}

/*!
 * Returns the current options of the solver.
 *
 * @return Solver options.
 */
template <typename CAESReal>
typename SolverContextImpl<CAESReal>::Options ECHMET_CC SolverContextImpl<CAESReal>::options() const noexcept
{
	return m_options;
}

/*!
 * Sets new solver options.
 *
 * @param[in] options New options.
 *
 * @retval RetCode::OK Success.
 * @retval RetCode::E_INVALID_ARGUMENT Nonsensical option value was passed as the argument.
 */
template <typename CAESReal>
RetCode ECHMET_CC SolverContextImpl<CAESReal>::setOptions(const Options options) noexcept
{
	/* Changing thread-safetiness is not allowed */
	if ((m_options & Options::DISABLE_THREAD_SAFETY) != (options & Options::DISABLE_THREAD_SAFETY))
		return RetCode::E_INVALID_ARGUMENT;

	m_options = options;

	return RetCode::OK;
}

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_SOLVERCONTEXTIMPL_HPP
