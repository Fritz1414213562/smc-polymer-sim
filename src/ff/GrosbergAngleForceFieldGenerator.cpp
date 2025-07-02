#include "GrosbergAngleForceFieldGenerator.hpp"

std::unique_ptr<OpenMM::Force> GrosbergAngleForceFieldGenerator::generate() const
{
	const std::string potential_formula = fmt::format(
		"{id}_k * (1 - cos(theta - {id}_v0))", fmt::arg("id", ffgen_id_));

	auto angle_ff = std::make_unique<OpenMM::CustomAngleForce>(potential_formula);
    angle_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
	angle_ff->addPerAngleParameter(fmt::format("{}_k",  ffgen_id_));
	angle_ff->addPerAngleParameter(fmt::format("{}_v0", ffgen_id_));

	for (std::size_t idx = 0; idx < indices_vec_.size(); ++idx)
	{
		const indices_type& triplet = indices_vec_[idx];
		angle_ff->addAngle(triplet[0], triplet[1], triplet[2], {ks_[idx], v0s_[idx]});
	}
	return angle_ff;
}
