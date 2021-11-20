#ifndef ECHMET_SYSCOMP_SYSCOMP_H
#define ECHMET_SYSCOMP_SYSCOMP_H

#include <echmetelems.h>
#include <echmetmodule.h>

namespace ECHMET {
namespace SysComp {

class AnalyticalConcentration;
class IonicConcentration;
class Constituent;
class ContainedLigandIonicForm;
class IonicForm;
class InComplexForm;
class InConstituent;
class InLigandForm;
class InLigandGroup;

typedef Vec<AnalyticalConcentration *> AnalyticalConcentrationVec;
typedef Vec<IonicConcentration *> IonicConcentrationVec;
typedef Vec<Constituent *> ConstituentVec;
typedef Vec<ContainedLigandIonicForm> ContainedLigandIonicFormVec;
typedef Vec<IonicForm *> IonicFormVec;
typedef Vec<InComplexForm> InCFVec;
typedef Vec<InConstituent> InConstituentVec;
typedef Vec<InLigandForm> InLFVec;
typedef Vec<InLigandGroup> InLGVec;

/* Input data structures. Used to construct the internal system composition representation */

/*!
 * Allowed types of constituent
 */
ECHMET_ST_ENUM(ConstituentType) {
	INVALID = 0,	/*!< Invalid or unknown type of constituent. */
	NUCLEUS = 1,	/*!< Constituent is a nucleus. */
	LIGAND = 2	/*!< Constituent is a ligand. */
	ENUM_FORCE_INT32_SIZE(SysCompConstituentType)
};

/*!
 * Describes a ligand form bound to a complex nucleus.
 */
class InLigandForm {
public:
	FixedString *ligandName;		/*!< Unique identifier of the ligand. */
	int32_t charge;				/*!< Electric charge of the ligand. */
	uint32_t maxCount;			/*!< Maximum number of ligands of this kind with the electric charge \p charge bound to a nucleus. */
	RealVec *pBs;				/*!< Vector of consecutive complexation stability constants.
						     First item corresponds to the form with one ligand bound,
						     second to two ligands etc. Number of entries in this vector
						     shall be the same as \p maxCount. */
	RealVec *mobilities;			/*!< Vector of electrophoretic mobilities of the complexed forms in <tt>(m\.m/V/s) \. 1e9</tt>.
						     First item corresponds to the form with one ligand bound,
						     second to two ligands etc. Number of entries in this vector
						     shall be the same as \p maxCount. */

};
IS_POD(InLigandForm)

/*!
 * Describes a group of ligands bound to a complex nucleus.
 */
class InLigandGroup {
public:
	InLFVec *ligands;			/*!< Vector of ligand forms that make up the group. */
};
IS_POD(InLigandGroup)

/*!
 * Describes a complex form of a complex nucleus.
 */
class InComplexForm {
public:
	int32_t nucleusCharge;			/*!< Electric charge of the complex nucleus. This charge must be within \p chargeLow and \p chargeHigh
						     values of the given nucleus. */
	InLGVec *ligandGroups;			/*!< Vector of ligand groups that describe the complexation relations of the given nucleus. */

};
IS_POD(InComplexForm)

/*!
 * Describes a chemical constituent contained in the system.
 */
class InConstituent {
public:
	ConstituentType ctype;			/*!< Role of the constituent in the complexation scheme. */
	FixedString *name;			/*!< Unique identifier of the constituent. */
	int32_t chargeLow;			/*!< Lowest electric charge the constituent can attain. */
	int32_t chargeHigh;			/*!< Highest electric charge the constituent can attain. */
	RealVec *pKas;				/*!< Vector of acidobazic pKa constants.
						     Number of items in the vector shall be <tt>chargeHigh - chargeLow</tt>. */
	RealVec *mobilities;			/*!< Vector of electrophoretic mobilities of the constituent in <tt>(V\.s/m/m) \. 1e9</tt>.
						     Number of items in the vector shall be <tt>chargeHigh - chargeLow + 1</tt>.*/

	InCFVec *complexForms;			/*!< If the constituent is a nucleus, this vector describes its complexation relations. */
	ECHMETReal viscosityCoefficient;	/*!< Amount of contribution of the constituent to overall viscosity of the solution.
						     The value's dimension assumes that total constituent concentration is expressed
						     in <tt>mmol/dm<sup>3</sup></tt>. */
};
IS_POD(InConstituent)

/*!
 * Memory pools to use with \p CalculatedProperties initialized with <tt>initializeCalculatedPropertiesWithPools</tt>
 */
class CalculatedPropertiesPools {
public:
	ECHMETReal *ionicConcentrations;
	ECHMETReal *ionicMobilities;
};
IS_POD(CalculatedPropertiesPools)

/*!
 * Ligand ionic form contained within a specific complexation relation.
 */
class ContainedLigandIonicForm {
public:
	const Constituent *ligand;		/*!< The ligand */
	int charge;				/*!< Electric charge of the ligand */
};
IS_POD(ContainedLigandIonicForm)

/*!
 * Type of component described by \p IonicForm.
 */
ECHMET_ST_ENUM(IonicFormType) {
	CONSTITUENT = 0,			/*!< Component is a regular constituent. */
	H = 1,					/*!< Component is a hydroxonium ion. */
	OH = 2					/*!< Component is a hydroxyle ion. */
	ENUM_FORCE_INT32_SIZE(SysCompIonicFormType)
};

/*!
 * Describes a concrete ionic form that is present in the chemical system.
 */
class IonicForm {
public:
	const Constituent *nucleus;				/*!< Complex nucleus the ionic form is base upon. */
	const FixedString *name;				/*!< Unique name of the ionic form.
								     This may be \p NULL if the ionic form corresponds to <tt>H+</tt> or <tt>OH-</tt>.
								     Such information can be retrieved from \p IonicConcentration \p ifType field. */
	int32_t nucleusCharge;					/*!< Electric charge of the nucleus */
	ECHMETReal limitMobility;				/*!< Limit electrophoretic mobility of the ionic form expressed in <tt>(m\.m/V/s) \. 1e9</tt>. */
	int32_t totalCharge;					/*!< Total charge of this form */

	/* Sensible only for complexes */
	const Constituent *ligand;				/*!< Ligand that makes up the complex.
								     Set only if the ionic form is a complex, \p NULL otherwise. */
	int32_t ligandCharge;					/*!< Electric charge of the ligand.
								     Set only if the ionic form is a complex, minimum allowed integer value otherwise */
	int32_t ligandCount;					/*!< Number of ligands of kind \p ligand bound to the ionic form.
								     Zero when ionic form is not a complex */
	IonicForm *ancestor;					/*!< Pointer to the ionic form the complex is based on.
								     \p NULL if the ionic form is not a complex */
	ECHMETReal pB;						/*!< Consecutive complexation stability constant of the complex.
								     Set to NaN if the ionic form is not a complex */
	ContainedLigandIonicFormVec *containedLigandIFs;	/*!< Vector of all ionic forms contained in the ionic form */

	size_t ionicConcentrationIndex;				/*!< Index in the respective array with the ionic concentration of the ionic form */
	size_t ionicMobilityIndex;				/*!< Index in the respective array with the ionic mobility of the ionic form */

	IonicFormType ifType;					/*!< Type of the ionic form */
};
IS_POD(IonicForm)

/*!
 * Describes a chemical constituent
 */
class Constituent {
public:
	ConstituentType ctype;					/*!< Role of the constituent in the complexation relations */
	FixedString *name;					/*!< Unique name of the constituent */
	int32_t chargeLow;					/*!< Lowest electric charge the constituent can attain */
	int32_t chargeHigh;					/*!< Highest electric charge the constituent can attain */
	IonicFormVec *ionicForms;				/*!< Vector of ionic forms of this constituent */
	RealVec *pKas;						/*!< Vector of pKa constants. Number of items in this vector shall be <tt>chargeHigh - chargeLow</tt>  */
	RealVec *limitMobilities;				/*!< Vector of limit electrophoretic mobilities in <tt>(m\.m/V/s) \. 1e9</tt>. Size of this vector shall be <tt>chargeHigh - chargeLow + 1</tt>. */

	size_t analyticalConcentrationIndex;			/*!< Index in the respective array with analytical concentration of the constituent */
	size_t effectiveMobilityIndex;				/*!< Index in the respective array with effective mobility of the constituent */

	ECHMETReal viscosityCoefficient;			/*!< Amount of contribution of the constituent to overall viscosity of the solution.
								     The value's dimension assumes that total constituent concentration is expressed
								     in <tt>mmol/dm<sup>3</sup></tt>. */
};
IS_POD(Constituent)

typedef SKMap<size_t> ConcentrationSKMap;

/*!
 * Description of a chemical system
 */
class ChemicalSystem {
public:
	ConstituentVec *constituents;				/*!< Vector of all constituents */
	IonicFormVec *ionicForms;				/*!< Vector of all ionic forms that can exist in the system */

	ConcentrationSKMap *analyticalConcentrationsByName;	/*!< Mapping of constituents' names to indices in the analytical concentrations array */
	ConcentrationSKMap *effectiveMobilitiesByName;		/*!< Mapping of constituents' names to indices in the effective mobilities array */
	ConcentrationSKMap *ionicConcentrationsByName;		/*!< Mapping of ionic forms' names to indices in the ionic concentrations array */
	ConcentrationSKMap *ionicMobilitiesByName;		/*!< Mapping of ionic forms' names to indices in the ionic mobilities array */
};
IS_POD(ChemicalSystem)

/*!
 * Properties of a given chemical system with given analytical concentrations
 * of its constituents.
 */
class CalculatedProperties {
public:
	RealVec *ionicConcentrations;				/*!< Vector of ionic concentrations sorted in order described by <tt>analyticalConcentrationsByName</tt> map */
	RealVec *ionicMobilities;				/*!< Vector of ionic mobilities sorted in order described by <tt>ionicMobilitiesByName</tt> map.
								     Note that these can be corrected for ionic strength using the Onsager-Fuoss law by a call to
								     <tt>correctMobilities()</tt> in the <tt>IonProps</tt>. */
	RealVec *effectiveMobilities;				/*!< Vector of ionic mobilities sorted in order described by <tt>effectiveMobilitiesByName</tt> map.
								     Note that these are calculated by a call to <tt>calculateEffectiveMobilities()</tt> in the <tt>IonProps</tt> package. */

	ECHMETReal ionicStrength;				/*!< Ionic strength of the system in <tt>mol/dm3</tt> */
	ECHMETReal conductivity;				/*!< Electric onductivity of the system in <tt>S/m</tt>.
								     Note that this is calculated by a call to <tt>calculateConducvitivy()</tt> in the <tt>IonProps</tt> package. */
};
IS_POD(CalculatedProperties)

extern "C" {

/*!
 * Initialize memory pools for blocks of <tt>ionicConcentrations</tt> and <tt>ionicMobilities</tt> array.
 *
 * The intended use of this function is to allocate memory space for multiple <tt>ionicConcentrations</tt> and <tt>ionicMobilities</tt>
 * vectors in a large contiguous memory space. The reason to do this is potentially better cache efficiency.
 * Memory allocated by the pool must be released by <tt>releaseCalculatedPropertiesPools()</tt>.
 *
 * @param[in,out] pools The \p CalculatedPropertiesPools object to be initialized
 * @param[in] chemSystem \p ChemicalSystem to initialize the \p CalculatedPropertiesPools for
 * @param[in] numBlocks Number of \p CalculatedProperties objects that are expected to be used simultaneously
 *
 * @retval OK Success
 * @retval E_NO_MEMORY Insufficient memory to allocate pools
 * @retval E_INVALID_ARGUMENT Number of blocks or ionic forms in the system is less than 1
 */
ECHMET_API RetCode ECHMET_CC allocateCalculatedPropertiesPools(CalculatedPropertiesPools &pools, const ChemicalSystem &chemSystem, const size_t numBlocks);

/*!
 * \brief Compares two \p InConstituent s
 *
 * Comparison is done on all members of the \p InConstituent class.
 * Two instances of this class are considered equal only if all
 * members of the class including their content are equal.
 *
 * If comparison of complexations is enabled for nuclei complexations
 * are considered equal only if they are laid out in the exactly same order
 * in both constituents.
 *
 * @param[in] first First \p InConstituent to compare
 * @param[in] second Second \p InConstituent to compare
 * @param[in] compareComplexations If true, even the complexation part of constituents is compared.
 *                                 This applies only to constituents of \p NUCLEUS type.
 *
 * @return True if the constituents are identical, false otherwise
 */
ECHMET_API bool ECHMET_CC compareInConstituents(const InConstituent &first, const InConstituent &second, const bool compareComplexations = true) ECHMET_NOEXCEPT;

/*!
 * \brief Compares two \p InLigandForm s
 *
 * Comparison is done on all members of the \p InLigandForm class.
 * Two instances of this class are considered equal only if all
 * members of the class including their content are equal.
 *
 * @param[in] first First \p InLigandForm to compare
 * @param[in] second Second \p InLigandForm to compare
 *
 * @return True if the forms are identical, false otherwise
 */
ECHMET_API bool ECHMET_CC compareInLigandForms(const InLigandForm &first, const InLigandForm &second) ECHMET_NOEXCEPT;

/*!
 * Returns an \p ECHMETVec of \p InComplexForm s
 *
 * @param[in] reserve Number of items to reserve memory for.
 *                    Zero forgoes the reservation
 * @return Pointer to vector on success, \p NULL on failure
 */
ECHMET_API InCFVec * ECHMET_CC createInCFVec(const size_t reserve) ECHMET_NOEXCEPT;

/*!
 * Returns an \p ECHMETVec of \p InConstituent s
 *
 * @param[in] reserve Number of items to reserve memory for.
 *                    Zero forgoes the reservation
 * @return Pointer to vector on success, \p NULL on failure
 */
ECHMET_API InConstituentVec * ECHMET_CC createInConstituentVec(const size_t reserve) ECHMET_NOEXCEPT;

/*!
 * Returns an \p ECHMETVec of \p InLigandForm s
 *
 * @param[in] reserve Number of items to reserve memory for.
 *                    Zero forgoes the reservation
 * @return Pointer to vector on success, \p NULL on failure
 */
ECHMET_API InLFVec * ECHMET_CC createInLFVec(const size_t reserve) ECHMET_NOEXCEPT;

/*!
 * Returns an \p ECHMETVec of \p InLigandGroup s
 *
 * @param[in] reserve Number of items to reserve memory for.
 *                    Zero forgoes the reservation
 * @return Pointer to vector on success, \p NULL on failure
 */
ECHMET_API InLGVec * ECHMET_CC createInLGVec(const size_t reserve) ECHMET_NOEXCEPT;

/*!
 * Deep-copies \p InConstituent object
 *
 * @param[out] dst \p InConstituent object
 * @param[in] src \p InConstituent object to be duplicated
 *
 * @retval OK Duplication successful
 * @retval E_NO_MEMORY No memory to duplicate \p InConstituent
 */
ECHMET_API RetCode ECHMET_CC duplicateInConstituent(InConstituent &dst, const InConstituent &src) ECHMET_NOEXCEPT;

/*!
 * Deep-copies \p InConsituentVec object
 *
 * @param[in] src \p InConstituentVec object to be copied
 *
 * @return Pointer to the duplicated vector, \p NULL on failure
 */
ECHMET_API InConstituentVec * ECHMET_CC duplicateInConstituentVec(const InConstituentVec *src) ECHMET_NOEXCEPT;

/*!
 * Initialize \p CalculatedProperties object
 *
 * @param[in,out] calcProps \p CalculatedPropertes to be initialized
 * @param[in] chemSystem \p ChemicalSystem to initialize the \p CalculatedProperties for
 */
ECHMET_API RetCode ECHMET_CC initializeCalculatedProperties(CalculatedProperties &calcProps, const ChemicalSystem &chemSystem) ECHMET_NOEXCEPT;

/*!
 * Initialize \p CalculatedProperties object using \p CalculatedPropertiesPools for <tt>ionicConcentrations</tt> and <tt>ionicMobilities</tt> vectors.
 * \p CalculatedProperties allocated by this function must be released with <tt>releaseCalculatedPropertiesWithPools</tt>.
 *
 * @param[in,out] calcProps \p CalculatedProperties to be initialized
 * @param[in] chemSystem \p ChemicalSystem to initialize the \p CalculatedProperties for
 * @param[in] pools \p CalculatedPropertiesPools to use
 * @param[in[ block Index of the \p CalculatedProperties object. Note that multiple \p CalculatedProperties objects used simultaneously must not use the same index
 */
ECHMET_API RetCode ECHMET_CC initializeCalculatedPropertiesWithPools(CalculatedProperties &calcProps, const ChemicalSystem &chemSystem, CalculatedPropertiesPools &pools, const int32_t block);

/*!
 * Create a vector with the correct size to hold analytical concentrations of all constituents
 * in a system
 *
 * @param[in,out] acVec Reference to the vector to be created
 * @param[in] chemSystem Chemical system to create the analytical concentrations vector for
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to create the vector
 */
ECHMET_API RetCode ECHMET_CC makeAnalyticalConcentrationsVec(RealVec *&acVec, const ChemicalSystem &chemSystem) ECHMET_NOEXCEPT;

/*!
 * Initializes \p ChemicalSystem and \p CalculatedProperties for a given input system composition.
 *
 * @param[in,out] chemSystem Reference to a \p ChemicalSystem object to be initialized
 * @param[in,out] calcProps Reference to a \p CalculatedProperties object to be initialized
 * @param[in] inputData Vector of input constituents that make up the chemical system
 *
 * @retval RetCode::OK Success
 * @retval RetCode::E_NO_MEMORY Insufficient memory to perform the initialization
 * @retval RetCode::E_INVALID_ARGUMENT Invalid object was passed as \p InConstituentVec
 * @retval RetCode::E_INVALID_CONSTITUENT System contains a constituent with invalid properties
 * @retval RetCode::E_DATA_TOO_LARGE System is too large to process
 * @retval RetCode::E_INVALID_COMPLEXATION System contains nonsensical complexation relation
 * @retval RetCode::E_DUPLICIT_CONSTITUENTS System contains multiple constituents with the same name
 */
ECHMET_API RetCode ECHMET_CC makeComposition(ChemicalSystem &chemSystem, CalculatedProperties &calcProps, const InConstituentVec *inputData) ECHMET_NOEXCEPT;

/*!
 * Convenience function that releases all resources claimed by
 * \p ChemicalSystem
 *
 * Make sure that no data from the given \p ChemicalSystem object are referenced
 * anywhere before calling this function.
 *
 * @param[in] chemSystem \p ChemicalSystem to be released
 */
ECHMET_API void ECHMET_CC releaseChemicalSystem(ChemicalSystem &chemSystem) ECHMET_NOEXCEPT;

/*!
 * Convenience function that releases all resources claimed by
 * \p CalculatedProperties initialized with <tt>initializeCalculatedProperties</tt> or <tt>makeComposition</tt>
 *
 * @param[in] calcProps \p CalculatedProperties to be released
 */
ECHMET_API void ECHMET_CC releaseCalculatedProperties(CalculatedProperties &calcProps) ECHMET_NOEXCEPT;

/*!
 * Frees memory allocated by <tt>allocateCalculatedPropertiesPools</tt>
 *
 * @param[in,out] pools \p CalculatedPropertiesPools object to be released
 */
ECHMET_API void ECHMET_CC releaseCalculatedPropertiesPools(CalculatedPropertiesPools &pools);

/*!
 * Convenience function that releases all resources claimed by
 * \p CalculatedProperties initialized with <tt>initializeCalculatedPropertiesWithPools</tt>
 *
 * @param[in] calcProps \p CalculatedProperties to be released
 */
ECHMET_API void ECHMET_CC releaseCalculatedPropertiesWithPools(CalculatedProperties &calcProps);

/*!
 * Convenience function to release \p InConstituent object
 *
 * @param[in] inC \p InConstituent object to be released
 */
ECHMET_API void ECHMET_CC releaseInConstituent(const InConstituent &inC) ECHMET_NOEXCEPT;

/*!
 * Convenience function that releases all resources claimed
 * by a vector of \p InConstituent s.
 *
 * @param[in] inVec Vector of \p InConstituent s to be disposed of.
 */
ECHMET_API void ECHMET_CC releaseInputData(const InConstituentVec *inVec) ECHMET_NOEXCEPT;

} // extern "C"

} // namespace SysComp
} // namespace ECHMET

#endif // ECHMET_SYSCOMP_SYSCOMP_H
