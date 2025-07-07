#pragma once
#ifndef NUCLEUS_NUCLEUSCONSTANTS_HPP_
#define NUCLEUS_NUCLEUSCONSTANTS_HPP_

namespace NucleusConstants {
constexpr const double nbTofm2 = 1.0e-7;
constexpr const double fm2Tonb = 1.0 / nbTofm2;

constexpr const double hbarc = 0.1973;  // GeV fm-1

constexpr const double GeVTofmm1 = 1.0 / hbarc;
constexpr const double fmToGeVm1 = 1.0 / hbarc;
constexpr const double GeVm1Tofm = hbarc;
constexpr const double fmm1ToGeV = hbarc;

// source: https://arxiv.org/pdf/2501.14872
constexpr const double lognorm_sigma = 0.637;  // unit 1
}  // namespace NucleusConstants

#endif  // NUCLEUS_NUCLEUSCONSTANTS_HPP_
