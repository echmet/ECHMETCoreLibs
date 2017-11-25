#ifndef ECHMET_CAES_TYPES_H
#define ECHMET_CAES_TYPES_H

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/MPRealSupport>

#define ECHMET_IMPORT_INTERNAL
#include <echmetelems.h>
#include <echmetsyscomp.h>
#undef ECHMET_IMPORT_INTERNAL


namespace ECHMET {
namespace CAES {

template <typename CAESReal>
using SolverMatrix =  Eigen::Matrix<CAESReal, Eigen::Dynamic, Eigen::Dynamic, static_cast<int>(Eigen::ColMajor) | static_cast<int>(Eigen::AutoAlign)>;
template <typename CAESReal>
using SolverVector = Eigen::Matrix<CAESReal, Eigen::Dynamic, 1, static_cast<int>(Eigen::ColMajor) | static_cast<int>(Eigen::AutoAlign)>;

class InvalidComplexationConfigurationException : public std::exception {
public:
	explicit InvalidComplexationConfigurationException() : std::exception() {}
	const char * what() const noexcept
	{
		return "Maximum number of bound ligands does not match the number of stability constants";
	}
};

template <typename CAESReal>
class Ligand;

/*!
 * Data structure that represents one ionic form of a ligand.
 */
template <typename CAESReal>
class LigandIonicForm {
public:
	explicit LigandIonicForm(const std::string &name, const int charge, const size_t index, const Ligand<CAESReal> *base);
	explicit LigandIonicForm(const LigandIonicForm &other, const uint32_t max, const std::vector<CAESReal> &pBs);

	bool operator==(const LigandIonicForm &other) const;

	const std::string name;			/*!< Ionic form name */
	const Ligand<CAESReal> *base;		/*!< Pointer to the ligand the ionic form is based on */
	const size_t index;			/*!< Global index in the vector of all ligand ionic forms */
	const int charge;			/*!< Electric charge of the ionic form */
	const uint32_t max;			/*!< Maximum number of the ligand atoms or molecules that can be bound to the analyte */
	const std::vector<CAESReal> pBs;	/*!< List of consecutive stability constants */
};

template <typename CAESReal>
using LigandIonicFormVec = std::vector<LigandIonicForm<CAESReal> *>;
template <typename CAESReal>
using LigandIonicFormVecVec = std::vector<LigandIonicFormVec<CAESReal>>;

/*!
 * Data structure that represents a ligand contained in the background
 * electrolyte. A ligand can form complexes with a complexing component.
 */
template <typename CAESReal>
class Ligand {
public:
	explicit Ligand(const std::string &name, const int chargeLow, const int chargeHigh,
			const size_t analyticalConcentrationIndex, const std::vector<CAESReal> &pKas);

	const std::string name;							/*!< Name of the ligand */
	const int chargeLow;							/*!< Lowest electric charge the ligand can have */
	const int chargeHigh;							/*!< Highest electric charge the ligand can have */
	const size_t analyticalConcentrationIndex;				/*!< Total concentration of the ligand in the system */
	const std::vector<CAESReal> pKas;					/*!< List of acidobazic dissociation constants */
	LigandIonicFormVec<CAESReal> ionicForms;				/*!< List of all ionic forms of the ligand */
};

template <typename CAESReal>
using LigandVec = std::vector<Ligand<CAESReal> *>;

template <typename CAESReal>
class Form;

/*!
 * Data structure representing a chemical species that can be contained either
 * in the background electrolyte or in the sample. The species can form single-core
 * complexes with an arbitrary number or type of ligands (some limitations apply,
 * see the algorithm description for more details).
 */
template <typename CAESReal>
class ComplexNucleus {
public:
	explicit ComplexNucleus(const std::string &name, const int n, const int p,
				const size_t analyticalConcentrationIndex,
				const std::vector<CAESReal> &pKas);
	~ComplexNucleus();

	const std::string name;							/*!< Name of the analyte */
	const int chargeLow;							/*!< Lowest electric charge the analyte can have */
	const int chargeHigh;							/*!< Highest electric charge the analyte can have */
	const size_t analyticalConcentrationIndex;				/*!< Total concentration of the complexing component */
	const std::vector<CAESReal> pKas;					/*!< List of acidobazic dissociation constants */
	std::vector<std::vector<Form<CAESReal> *>> forms;			/*!< List of all complex forms the complexing component forms at all its ionic states */
};

template <typename CAESReal>
using CNVec = std::vector<ComplexNucleus<CAESReal> *>;

/*!
 * Data structure representing a ligand ionic form contained in a particular complex form.
 */
template <typename CAESReal>
class ContainedLigandIonicForm {
public:
	explicit ContainedLigandIonicForm();
	explicit ContainedLigandIonicForm(const uint32_t count, const LigandIonicForm<CAESReal> *lIF);

	ContainedLigandIonicForm &operator=(const ContainedLigandIonicForm &other);

	uint32_t count;					/*!< Number of ligands of the given ionic form bound in the given complex form */
	const LigandIonicForm<CAESReal> * const lIF;	/*!< Pointer to the ligand ionic form */
};

template <typename CAESReal>
using ContainedLigandIonicFormVec = std::vector<ContainedLigandIonicForm<CAESReal>>;

/*!
 * Data structure that represents a complex form of one complexing compoment with
 * ligands.
 */
template <typename CAESReal>
class Form {
public:
	explicit Form(const std::string &name, const int totalCharge,
		      const size_t ligandIFIdx = SIZE_MAX, const ContainedLigandIonicFormVec<CAESReal> &ligandsContained = ContainedLigandIonicFormVec<CAESReal>(),
		      const size_t ancestorIdx = SIZE_MAX, const size_t ancestorGlobalIdx = SIZE_MAX,
		      const size_t myIdx = 0, const CAESReal pB = nan(""));

	const std::string name;					/*!< Name of the complex form */
	const int totalCharge;					/*!< Total charge of the form */
	const size_t ligandIFIdx;				/*!< Index of the ligand ionic form that makes up this from in the vector of all ligand ionic forms.
								 * If a form contains more ligands, then this is the ligand that binds to the ancestor form like this:
								 * <tt>ML(0) + X(0) --> ML(0)X(0)</tt>, this value would then contain the index of ligand ionic form <tt>X(0)</tt>.
								 * Value of SIZE_MAX indicates a form that contains no ligand. */
	const size_t myIdx;					/*!< Index of the form in the vector of all forms belonging to a given complexing component */
	const size_t ancestorIdx;				/*!< Index of the ancestor form. Value of <tt>SIZE_MAX</tt> indicates a form that has no ancestor */
	const size_t ancestorGlobalIdx;				/*!< Index of the ancestor form. Value of <tt>SIZE_MAX</tt> indicates a form that has no ancestor */
	const ComplexNucleus<CAESReal> *cn;			/*!< Pointer to the complex nucleus the form is built upon */
	const ContainedLigandIonicFormVec<CAESReal> ligandsContained;	/*!< Vector of all ligand ionic forms contained in this form */
	CAESReal pB;						/*!< Value of pB for this complex form */
};

template <typename CAESReal>
using FormVec =  std::vector<Form<CAESReal> *>;

} // namespace CAES
} // namespace ECHMET

#include "types.hpp"

#endif // ECHMET_CAES_TYPES_H

