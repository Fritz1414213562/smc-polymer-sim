#ifndef POLYMER_MDMC_CONSTANTS_HPP
#define POLYMER_MDMC_CONSTANTS_HPP

#include <OpenMM.h>

#include <array>
#include <map>
#include <string>
#include <utility>

#include <cmath>

namespace Constant
{
// define physical constnt
inline static const double Na   = 6.02214076e23;         // avogadro constant [/mol]
inline static const double kB   = 1.380649e-23;          // boltzmann constant [J/K]
inline static const double pi   = 3.1415926535897932385;
inline static const double eps0 = 8.854187817e-12 /*[F/m] == [C^2/J/m]*/ * 1.0e3 /*[/J]->[/KJ]*/ * 1.0e-9 /*[/m]->[/nm]*/ / Na; // vacuum permittivity [C^2 mol/KJ/nm]
inline static const double elementary_charge = 1.6021766208e-19; // [C]
// define forcefield specific parameters
//constexpr double prexv_rmin12 = std::sqrt(6.0 / 7.0);
constexpr double prexv_rmin12 = 0.9258200997725514;// = std::sqrt(6.0 / 7.0);
constexpr double prexv_emin12 = 46656.0 / 823543.0;

}


#endif // POLYMER_MDMC_CONSTANTS_HPP
