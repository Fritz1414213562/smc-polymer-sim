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
			const std::vector<double> bond_lengths, const std::vector<double> sigmas,
			const std::vector<double> phi0s, const std::vector<double> theta0s,
			const bool use_periodic);

		std::unique_ptr<OpenMM::Force> generate() const override;

		const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }

		std::string name() const override { return "SegmentParallelization"; }
	
	private:
		std::vector<indices_type> indices_vec_;
		std::vector<double>       bond_ks_;
		std::vector<double>       dihedral_ks_;
		std::vector<double>       bond_lengths_;
		std::vector<double>       sigmas_;
		std::vector<double>       phi0s_;
		std::vector<double>       theta0s_;
		bool                      use_periodic_;
		std::string               ffgen_id_;
		std::string potential_formula_ =
			"- {id}_bond_k * f_r * g_theta * g_phi;"
			"g_theta = 1.0 - {id}_dihd_k * sin(theta - {id}_theta0)^2;"
			"g_phi   = 1.0 - {id}_dihd_k * sin(phi   - {id}_phi0)^2;"
 			//"f_r       = exp(-dr^2/ (2 * {id}_sigma * {id}_sigma));"
 			"f_r       = -1 / ({id}_sigma * {id}_sigma) * dr^2;"
			"dr = distance(p1, p3) - {id}_r0;"
			"phi = angle(p2, p1, p3) - angle(p1, p3, p4);"
			"theta = dihedral(p2, p1, p3, p4);"
			"pi = 3.1415926535897932385;";
 			//"g_theta   = max(g1*rect0t, rect1t);"
 			//"g_phi     = max(g2*rect0p, rect1p);"
 			//"rect0t    = step(pi   + dtheta) * step(pi   - dtheta);"
 			//"rect1t    = step(pi/2 + dtheta) * step(pi/2 - dtheta);"
 			//"rect0p    = step(pi   + dphi)   * step(pi   - dphi);"  
 			//"rect1p    = step(pi/2 + dphi)   * step(pi/2 - dphi);"  
			//"g1 = 1 - cos(dtheta)^2;"
			//"g2 = 1 - cos(dphi)^2;"
			//"dtheta = K * (theta - {id}_theta0);"
			//"dphi = K * (phi - {id}_phi0);"
			//"K = pi/(2 * {id}_dihd_k);"
			//"dr = distance(p1, p3) - {id}_r0;"
			//"phi = angle(p2, p1, p3) - angle(p1, p3, p4);"
			//"theta = dihedral(p2, p1, p3, p4);"
			//"pi = 3.1415926535897932385;";

			//"Ub + Ua + Ud;"
			//"Ud = {id}_dihd_k * (1 - cos(2 * (theta - {id}_theta0)));"
			//"Ua = {id}_dihd_k * (1 - cos(2 * (phi - psi - {id}_phi0)));"
			//"Ub = {id}_bond_k * (r - {id}_r0)^2;"
			//"Ub = {id}_bond_k * (r - {id}_r0)^2 * rect_r;"
			//"rect_r = step(r - {id}_r0);" // 1 when r > r0, else 0
			//"phi = angle(p2, p1, p3);"
			//"psi = angle(p1, p3, p4);"
			//"r = distance(p1, p3);"
			//"theta = dihedral(p2, p1, p3, p4);";

};


#endif // POLYMER_MDMC_SEGMENT_PARALLELIZATION_FORCE_FIELD_GENERATOR_HPP
