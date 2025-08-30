#include "GrosbergRepulsiveForceFieldGenerator.hpp"
#include "InteractionGroup.hpp"
#include "ForceFieldIDGenerator.hpp"
#include "src/util/Constants.hpp"

#include <iostream>
#include <limits>

// The size of the vector representing the per-particle parameters
// (in this case, radii) must match the system size. This must be guaranteed
// inside the read function.
GrosbergRepulsiveForceFieldGenerator::GrosbergRepulsiveForceFieldGenerator(
	const double eps,
	const std::vector<std::optional<double>>& radii,
	const index_pairs_type& ignore_list, const bool use_periodic,
	const std::vector<std::pair<std::string, std::string>> ignore_group_pairs,
	const std::vector<std::optional<std::string>> group_vec)
	: eps_(eps), radii_(radii),
	  ignore_list_(ignore_list), use_periodic_(use_periodic),
	  ffgen_id_(fmt::format("GBEXV{}", ffid.gen()))
{
	interaction_groups_ = extract_interaction_group(radii_, ignore_group_pairs, group_vec);
}

std::unique_ptr<OpenMM::Force> GrosbergRepulsiveForceFieldGenerator::generate() const
{
	const std::string potential_formula = fmt::format(
		"U * rect_r;"
		"rect_r = step(sigma_sum * {id}_sigma_factor - r);"
		"U = {id}_epsilon * (r12 - r6) + {id}_epsilon;"
		"r12 = r6 * r6;"
		"r6 = rrev * rrev * rrev * rrev * rrev * rrev;"
		"rrev = sigma_sum / r;"
		"sigma_sum = {id}_sigma1 + {id}_sigma2;",
		fmt::arg("id", ffgen_id_));

	auto exv_ff = std::make_unique<OpenMM::CustomNonbondedForce>(potential_formula);

	exv_ff->addPerParticleParameter(fmt::format("{}_sigma", ffgen_id_));
	exv_ff->addGlobalParameter(fmt::format("{}_epsilon", ffgen_id_), eps_);
	exv_ff->addGlobalParameter(fmt::format("{}_sigma_factor", ffgen_id_), std::pow(2.0, 1.0/6.0));

	double max_radius = std::numeric_limits<double>::min();
	double second_max_radius = std::numeric_limits<double>::min();
	for (std::size_t idx = 0; idx < radii_.size(); ++idx)
	{
		const std::optional<double>& radius = radii_[idx];
		if (radius)
		{
			const double radius_val = radius.value();
			exv_ff->addParticle({radius_val});
			if (max_radius < radius_val)
			{
				second_max_radius = max_radius;
				max_radius = radius_val;
			}
			else if (second_max_radius < radius_val)
			{
				second_max_radius = radius_val;
			}
		}
		else
		{
			exv_ff->addParticle({std::numeric_limits<double>::quiet_NaN()});
		}
	}

	// if interaction_groups size is 0, no interaction group will be added,
	// so all the particle in the system will be considerd as participant
	for(const auto& group_pair : interaction_groups_)
	{
		exv_ff->addInteractionGroup(group_pair.first, group_pair.second);
	}

	// set pbc condition
	if(use_periodic_)
	{
		exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
	}
	else
	{
		exv_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffNonPeriodic);
	}

	// set cutoff
	const double cutoff_distance   = (max_radius + second_max_radius)*std::pow(2.0, 1.0/6.0);
	std::cerr << "	  GrosbergRepulsiveExcludedVolume      : cutoff disntace is "
			  << cutoff_distance << " nm" << std::endl;
	exv_ff->setCutoffDistance(cutoff_distance);

	//const double cutoff_correction = std::pow(1.0 / cutoff_distance, 12);
	//exv_ff->addGlobalParameter(fmt::format("{}_cutoff_correction", ffgen_id_), cutoff_correction);

	// set exclusion list
	for(const auto& pair : ignore_list_)
	{
		exv_ff->addExclusion(pair.first, pair.second);
	}

	return exv_ff;
}


