#ifndef POLYMER_MDMC_HARMONIC_BOND_FORCE_FIELD_GENERATOR_HPP
#define POLYMER_MDMC_HARMONIC_BOND_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include <OpenMM.h>
#include <fmt/core.h>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

// [Caution!] This potential formulaw is "1/2*k*(r - r0)^2", so interaction
// coefficient is devided by 2. So if you use the coefficient base on
// "k*(r  - r0)^2" potential, you should double that value.
class HarmonicBondForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    HarmonicBondForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& v0s, const std::vector<double>& ks,
        const bool use_periodic)
        : indices_vec_(indices_vec), v0s_(v0s), ks_(ks), use_periodic_(use_periodic)
    {
        if(!(indices_vec.size() == v0s.size() && v0s.size() == ks.size()))
        {
            std::ostringstream oss;
            oss << "[error] HarmonicBondForceFieldGenerator: "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "v0 ("          << v0s.size()         << ") and "
                   "k ("           << ks.size()          << ") is not matched."
                << "The number os these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override;

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string name() const override { return "HarmonicBond"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       v0s_;
    std::vector<double>       ks_;
    bool                      use_periodic_;
};

#endif // POLYMER_MDMC_HARMONIC_BOND_FORCE_FIELD_GENERATOR_HPP
