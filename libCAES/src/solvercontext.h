#ifndef SOLVERCONTEXT_H
#define SOLVERCONTEXT_H

#include "types.h"

namespace ECHMET {
namespace CAES {

class SolverContext {
public:
	explicit SolverContext(const LigandVec *allLigands, const LigandIonicFormVec *allLigandIFs,
			       const CNVec *complexNuclei, const FormVec *allForms,
			       const SolverMatrix *preJacobian, const size_t concentrationCount) noexcept;
	~SolverContext() noexcept;

	const LigandVec *allLigands;			/*!< Vector of all ligands */
	const LigandIonicFormVec *allLigandIFs;		/*!< Vector of all ligand ionic forms */
	const CNVec *complexNuclei;			/*!< Vector of all complexing component */
	const FormVec *allForms;			/*!< Vector of all generated complex forms */
	const arma::mat *preJacobian;			/*!< Pregenerated part of the Jacobian */
	const size_t concentrationCount;		/*!< Total number of concentrations to be calculated */

};

} // namespace CAES
} // namespace ECHMET

#endif // SOLVERCONTEXT_H
