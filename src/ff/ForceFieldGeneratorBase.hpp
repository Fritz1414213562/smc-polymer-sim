#ifndef POLYMER_MDMC_FORCE_FIELD_GENERATOR_BASE_HPP
#define POLYMER_MDMC_FORCE_FIELD_GENERATOR_BASE_HPP

#include <OpenMM.h>

#include <memory>
#include <string>

class ForceFieldGeneratorBase
{
  public:
	virtual ~ForceFieldGeneratorBase() = default;
    virtual std::unique_ptr<OpenMM::Force> generate() const = 0;
    virtual std::string                    name()     const = 0;
};

#endif // POLYMER_MDMC_FORCE_FIELD_GENERATOR_BASE_HPP
