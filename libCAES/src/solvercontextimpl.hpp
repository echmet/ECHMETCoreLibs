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
					       const size_t analyticalConcentrationCount) noexcept :
	allLigands(allLigands),
	allLigandIFs(allLigandIFs),
	complexNuclei(complexNuclei),
	allForms(allForms),
	preJacobian(preJacobian),
	concentrationCount(concentrationCount),
	analyticalConcentrationCount(analyticalConcentrationCount)
{
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

} // namespace CAES
} // namespace ECHMET

#endif // ECHMET_CAES_SOLVERCONTEXTIMPL_HPP
