#include "ionprops_p.h"

namespace ECHMET {
namespace IonProps {

ECHMETReal ECHMET_CC calculateConductivity(const ComputationContext *ctx, SysComp::CalculatedProperties &calcProps) noexcept
{
	return calculateConductivityInternal<ECHMETReal>(ctx, calcProps);
}

RetCode ECHMET_CC calculateEffectiveMobilities(const ComputationContext *ctx, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	return calculateEffectiveMobilitiesInternal<ECHMETReal>(ctx, analyticalConcentrations, calcProps);
}

ECHMETReal ECHMET_CC calculatepH(const ComputationContext *ctx, const NonidealityCorrections corrections, SysComp::CalculatedProperties &calcProps) noexcept
{
	return calculatepHInternal<ECHMETReal>(ctx, corrections, calcProps);
}

ECHMETReal ECHMET_CC calculatepH_direct(const ECHMETReal &cH, const ECHMETReal &ionicStrength) noexcept
{
	return calculatepH_directInternal<ECHMETReal>(cH, ionicStrength);
}

RetCode ECHMET_CC correctMobilities(ComputationContext *ctx, const NonidealityCorrections corrections, const RealVec *analyticalConcentrations, SysComp::CalculatedProperties &calcProps) noexcept
{
	return correctMobilitiesInternal<ECHMETReal>(ctx, corrections, analyticalConcentrations, calcProps);
}

RetCode ECHMET_CC getTransferMultiplier(const SysComp::Constituent *c, const SysComp::IonicForm *iF, int &xfrMult) noexcept
{
	switch (c->ctype) {
	case SysComp::ConstituentType::NUCLEUS:
		xfrMult = 1;
		break;
	case SysComp::ConstituentType::LIGAND:
	{
		auto getLigandCount = [](const SysComp::IonicForm *iF, const FixedString *name) {
			for (;;) {
				if (iF->ligand == nullptr)
					return 1; /* This ionic form corresponds to free ligand */

				if (*(iF->ligand->name) == *name)
					return iF->ligandCount;

				if (iF->ancestor == nullptr)
					return 0;

				iF = iF->ancestor;
			}
		};

		xfrMult = getLigandCount(iF, c->name);
		if (xfrMult == 0)
			return RetCode::E_BAD_INPUT;
		}
		break;
	default:
		return RetCode::E_BAD_INPUT;
	}

	return RetCode::OK;
}

ComputationContext * ECHMET_CC makeComputationContext(const SysComp::ChemicalSystem &chemSystem, const ComputationContext::Options options) noexcept
{
	return new (std::nothrow) ComputationContextImpl<ECHMETReal>(chemSystem, options);
}

ComputationContext::~ComputationContext() noexcept {}

} // namespace IonProps

} // namespace ECHMET
