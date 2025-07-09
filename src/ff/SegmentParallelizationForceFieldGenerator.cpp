#include "SegmentParallelizationForceFieldGenerator.hpp"


SegmentParallelizationForceFieldGenerator::SegmentParallelizationForceFieldGenerator(
	const std::vector<indices_type>& indices_vec, const std::vector<double> bond_ks,
	const std::vector<double> dihedral_ks, const std::vector<double> bond_lengths,
	const std::vector<double> phases, const bool use_periodic)
	: indices_vec_(indices_vec), bond_ks_(bond_ks), dihedral_ks_(dihedral_ks),
	  bond_lengths_(bond_lengths), phases_(phases), use_periodic_(use_periodic),
	  ffgen_id_(fmt::format("SGP{}", ffid.gen()))
{
	if ( !(indices_vec_.size() == bond_ks_.size()
		&& indices_vec_.size() == dihedral_ks_.size()
		&& indices_vec_.size() == bond_lengths_.size()
		&& indices_vec_.size() == phases_.size()))
	{
		std::ostringstream oss;
		oss << "[error] SegmentParallelizationForceFieldGenerator: "
			   "parameter number of "
			   "indices_vec  (" << indices_vec_.size()  << "), "
			   "bond_ks  ("     << bond_ks_.size()      << "), "
			   "dihedral_ks  (" << dihedral_ks_.size()  << "), "
			   "bond_lengths (" << bond_lengths_.size() << "), and "
			   "phases       (" << phases_.size()       << ") is not matched."
			<< "The number of these parameters must be the same.";
		throw std::runtime_error(oss.str());
	}
}

std::unique_ptr<OpenMM::Force> SegmentParallelizationForceFieldGenerator::generate() const
{
	const std::string potential_formula = fmt::format(potential_formula_, fmt::arg("id", ffgen_id_));
	auto bond_ff = std::make_unique<OpenMM::CustomCompoundBondForce>(4, potential_formula);
	bond_ff->setUsesPeriodicBoundaryConditions(use_periodic_);
	bond_ff->addPerBondParameter(fmt::format("{}_bond_k", ffgen_id_));
	bond_ff->addPerBondParameter(fmt::format("{}_r0",     ffgen_id_));
	bond_ff->addPerBondParameter(fmt::format("{}_dihd_k", ffgen_id_));
	bond_ff->addPerBondParameter(fmt::format("{}_theta0", ffgen_id_));

	for (std::size_t idx = 0; idx < indices_vec_.size(); ++idx)
	{
		const std::array<std::size_t, 4>& pairs = indices_vec_[idx];
		bond_ff->addBond({pairs.begin(), pairs.end()},
			{bond_ks_[idx], bond_lengths_[idx], dihedral_ks_[idx], phases_[idx]});
	}
	return bond_ff;
}

