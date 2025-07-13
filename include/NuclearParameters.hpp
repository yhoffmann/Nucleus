/*
 * sources:
 *
 * [1]
 *     Improved version of the PHOBOS Glauber Monte Carlo
 *     Loizides, C. et al. SoftwareX, Volume 1, 13 - 18, 2015
 *
 * [2]
 *     Improved Monte Carlo Glauber predictions at present and future nuclear
 *     colliders Constantin Loizides, Jason Kamin, and David d'Enterria Phys.
 *     Rev. C 97, 054910, 2018
 */

#pragma once
#ifndef NUCLEUS_NUCLEARPARAMETERS_HPP_
#define NUCLEUS_NUCLEARPARAMETERS_HPP_

#include <stdlib.h>

#include <iostream>

#define NUCLEUS_ERROR_NO_NUC_PARAMS 30

namespace NuclearParameters {
static constexpr const char* const ERROR_NO_DATA =
    "No data for this nucleon number. Exiting!";

struct NucleusParameters {
  double mean_bulk_radius;            // R
  double mean_surface_diffusiveness;  // a
  double deviation_from_spherical;    // w
};

enum Type : uint {
  p = 1,
  O16 = 16,
  Si28 = 28,
  S32 = 32,
  Ca40 = 40,
  Ni58 = 58,
  Cu62 = 62,
  Cu63 = 63,
  Xe129 = 129,
  W186 = 186,
  Au197 = 197,
  Pb207 = 207,
  Pb208 = 208,
};

// data from [1] and [2]
const struct {
  NucleusParameters p = {1.0, 0.0001, 0.0};
  NucleusParameters O16 = {2.608, 0.513, -0.51};
  NucleusParameters Si28 = {3.34, 0.58, -0.233};
  NucleusParameters S32 = {2.54, 2.191, 0.16};
  NucleusParameters Ca40 = {3.766, 0.586, -0.161};
  NucleusParameters Ni58 = {4.309, 0.517, -0.1308};
  NucleusParameters Cu62 = {4.2, 0.596, 0.0};
  NucleusParameters Cu63 = Cu62;
  NucleusParameters Xe129 = {5.36, 0.59, 0.0};
  NucleusParameters W186 = {6.58, 0.48, 0.0};
  NucleusParameters Au197 = {6.38, 0.535, 0.0};
  NucleusParameters Pb207 = {6.62, 0.546, 0.0};
  NucleusParameters Pb208 = Pb207;
} params;

inline const NucleusParameters& get(uint atomic_num) {
  switch (atomic_num) {
    case p:
      return params.p;

    case O16:
      return params.O16;

    case Si28:
      return params.Si28;

    case S32:
      return params.S32;

    case Ca40:
      return params.Ca40;

    case Ni58:
      return params.Ni58;

    case Cu62:
      return params.Cu62;

    case Cu63:
      return params.Cu63;

    case Xe129:
      return params.Xe129;

    case W186:
      return params.W186;

    case Au197:
      return params.Au197;

    case Pb207:
      return params.Pb207;

    case Pb208:
      return params.Pb208;

    default:
      std::cerr << ERROR_NO_DATA << std::endl;
      exit(NUCLEUS_ERROR_NO_NUC_PARAMS);
  }
}
}  // namespace NuclearParameters

#endif  // NUCLEUS_NUCLEARPARAMETERS_HPP_
