#ifndef POLYMER_MDMC_SEGMENT_PARALLELIZATION_FORCE_FIELD_GENERATOR_HPP
#define POLYMER_MDMC_SEGMENT_PARALLELIZATION_FORCE_FIELD_GENERATOR_HPP
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"
#include <OpenMM.h>
#include <fmt/core.h>
#include <memory>
#include <sstream>
#include <string>
#include <array>

class SegmentParallelizationForceFieldGenerator final : public ForceFieldGeneratorBase
{
	public:
		using indices_type = std::array<std::size_t, 4>;
	public:
		SegmentParallelizationForceFieldGenerator(
			const std::vector<indices_type>& indices_vec,
			const std::vector<double> bond_ks, const std::vector<double> dihedral_ks,
			const std::vector<double> bond_lengths, const std::vector<double> phi_0s,
			const std::vector<double> psi_0s, const std::vector<double> phases,
			const bool use_periodic);

		std::unique_ptr<OpenMM::Force> generate() const override;

		const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

		std::string name() const override { return "SegmentParallelization"; }
	
	private:
		std::vector<indices_type> indices_vec_;
		std::vector<double>       bond_ks_;
		std::vector<double>       dihedral_ks_;
		std::vector<double>       bond_lengths_;
		std::vector<double>       phi_0s_;
		std::vector<double>       psi_0s_;
		std::vector<double>       phases_;
		bool                      use_periodic_;
		std::string               ffgen_id_;
		std::string potential_formula_ =
			"{id}_bond_k * (r - {id}_r0)^2 + {id}_bond_k * (phi - {id}_phi0)^2 + {id}_bond_k * (psi - {id}_psi0)^2 + {id}_dihd_k * (1 - cos(2 * (theta - {id}_theta0)));"
			"phi = angle(p2, p1, p3);"
			"psi = angle(p1, p3, p4);"
			"r = distance(p1, p3);"
			"theta = dihedral(p2, p1, p3, p4);";

};


#endif // POLYMER_MDMC_SEGMENT_PARALLELIZATION_FORCE_FIELD_GENERATOR_HPP
