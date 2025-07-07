#pragma once
#ifndef NUCLEUS_NUCLEUS_HPP_
#define NUCLEUS_NUCLEUS_HPP_

#include <math.h>
#include <stdlib.h>

#include <random>
#include <string>

#include "NucleusConstants.hpp"

// specifies which distribution to use for nucleon position sampling in the
// rejection sampling
enum SamplingDistribution : unsigned char { WoodsSaxon, Gaussian };

// contains 3d position of nucleon (basically just a double[3])
struct NucleonPos {
  double x, y, z;

  friend std::ostream& operator<<(std::ostream& stream, const NucleonPos& pos);
};

/* Nucleus class
 * -allows for the modelling a simple atomic nucleus using, for example, the
 * MC-Glauber model
 * -quick access to nucleon positions for collision/scattering simulations
 *
 * Construct using the given constructor. Access to nucleon positions is
 * possible with the method `get_nucleon_pos(uint)`.
 */
class Nucleus {
 protected:
  // atomic number
  uint m_atomic_num = 1;
  // 1d size parameter of nucleon, for example stddev of normal distribution
  double m_nucleon_size = std::sqrt(3.3);  // GeVm1
  // avg radius of nuclei
  double m_mean_bulk_radius;  // GeVm1
  // nucleus surface diffusiveness
  double m_mean_surface_diffusiveness;  // GeVm1
  // 3D positions of nucleons, relative to center of mass
  NucleonPos* m_nucleon_pos = nullptr;
  // weights for invidual nucleon thicknesses; is not automatically
  // sampled/initialized (use sample_nucleon_weights())
  double* m_nucleon_weights = nullptr;

  /* parameter describing range of rejection sampling values:
   * larger values are more accurate (with quickly deminishing returns),
   * smaller values result in faster sampling;
   * a good value is R + 10 * a
   */
  double m_sampling_range;
  // describes the nucleus model with simple rejection sampling to assert
  // nucleon postions
  SamplingDistribution m_sampling_distribution;
  // the random-number generator
  std::mt19937* m_rng = nullptr;
  // uniform real dist from 0 to 1, used for rejection sampling
  std::uniform_real_distribution<double> m_dist_01u =
      std::uniform_real_distribution<double>(0.0, 1.0);
  // lognorm dist for nucleon/hotspot weight sampling
  std::lognormal_distribution<double> m_dist_lognorm =
      std::lognormal_distribution<double>(0.0, NucleusConstants::lognorm_sigma);

  // the actual function returning numbers on U(0,1)
  inline double m_rand() { return m_dist_01u(*m_rng); }
  // inverse of the divisor needed to retain an average of 1 for the lognorm
  // dist
  const double lognorm_rand_inverse_divisor =
      1.0 / exp(NucleusConstants::lognorm_sigma *
                NucleusConstants::lognorm_sigma / 2.0);
  // the actual function returning lognormally distributed numbers with avg 1
  inline double m_lognorm_rand() {
    return m_dist_lognorm(*m_rng) * lognorm_rand_inverse_divisor;
  }

 public:
  // sample the nucleus or derived class (not including nucleon/hotspot weights)
  virtual void sample();
  // samples nucleon weights
  void sample_nucleon_weights();
  // sample nucleon weights (but fix nucleus total thicknes to always be the
  // same)
  void sample_nucleon_weights_fixed_avg();
  // get nucleon weight for nucleon n
  double get_nucleon_weight(uint n) const;
  // write nucleon positions to file
  void export_nucleon_positions(double impact_param_x, double impact_param_y,
                                const std::string& filepath) const;
  /* calculate and return nucleus thickness (nucleon density) at
   * transverse position (x,y)
   * NOTE: integrating this over the whole plane will yield the number of
   * nucleons, not 1
   */
  double get_nucleon_thickness(double x, double y) const;
  uint get_atomic_num() const;
  // get pointer to NucleonPos object containing pos of nucleon n
  const NucleonPos* get_nucleon_pos(uint nucleon_num) const;
  // memory-safe version of get_nucleon_pos (with check, thus slower)
  const NucleonPos* safe_get_nucleon_pos(uint nucleon_num) const;
  double get_mean_bulk_radius() const;
  double get_mean_surface_diffusiveness() const;
  // set size parameter of nucleon
  virtual void set_nucleon_size(double sigma_nn);
  double get_nucleon_size() const;
  // seed the rng with a uint
  void seed(uint seed);
  // copy-construct rng from another rng's current state
  void seed(const std::mt19937& rng);

  // removing default constructor, as we need the number of nucleons in order to
  // allocate the correct amount of memory
  Nucleus() = delete;
  // constructor will allocate memory and set all parameters based on the
  // specified nucleon number, and the sample the nucleon positions
  Nucleus(uint seed, uint atomic_num = 1, double nucleon_size = std::sqrt(3.3),
          SamplingDistribution sampling_distribution = WoodsSaxon);
  // copy constructor
  Nucleus(const Nucleus&);
  // move constructor
  Nucleus(Nucleus&&);
  // copy assignment
  Nucleus& operator=(const Nucleus&);
  // move assignment
  Nucleus& operator=(Nucleus&&);
  // virtual destructor
  virtual ~Nucleus();

 protected:
  void set_mean_bulk_radius();
  void set_mean_surface_diffusiveness();
  void set_sampling_range();
  void safe_delete_pos();
  void prepare_pos();
  void safe_delete_nucleon_weights();
  void prepare_nucleon_weights();
  void safe_delete_rng();
  void prepare_rng(uint seed);
  void prepare_rng(const std::mt19937& rng);
  void sample_single_pos(NucleonPos* nucleon_pos);
  bool fits_nucleon_distribution(double r_sqr);
};

/* Namespace containing sampling distributions (= nuclear models);
 * does currently not support correlations between positions;
 * does currently not support non-sphericality
 */
namespace NormedSamplingDistributions {
inline double woods_saxon(double r, double mean_bulk_radius,
                          double mean_surface_diffusiveness) {
  return (1.0 + exp(-mean_bulk_radius / mean_surface_diffusiveness)) /
         (1.0 + exp((r - mean_bulk_radius) / mean_surface_diffusiveness));
}

inline double gaussian(double r_sqr, double sigma_sqr) {
  return exp(-r_sqr / (2.0 * sigma_sqr));
}
}  // namespace NormedSamplingDistributions

#endif  // NUCLEUS_NUCLEUS_HPP_
