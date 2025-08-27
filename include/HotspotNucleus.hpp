#pragma once
#ifndef NUCLEUS_HOTSPOTNUCLEUS_HPP_
#define NUCLEUS_HOTSPOTNUCLEUS_HPP_

#include <iostream>
#include <random>

#include "Nucleus.hpp"

// contains transverse (!) position of hotspot relative to the center of mass of
// the whole nucleus
struct HotspotPos {
  double x, y;

  friend std::ostream& operator<<(std::ostream& stream, const HotspotPos& pos);
};

/* HotspotNucleus class
 * -derived class on top of a Nucleus object
 * -contains additionally the transverse positions of hotspots
 *
 * Construct this as you would the Nucleus class.
 * Access the hotspot positions with either `get_hotspot_pos(uint)` or
 * `get_hotspot_pos(uint, uint)`, where the former takes the absolute hotspot
 * index (e.g. `5` for the "last" hotspot in nucleon 1 for
 * `num_hotspots_per_nucleon=3`), and the latter takes both nucleon index and
 * hotspot index (using the same example as before: `1` and `2`).
 *
 * Also has support for giving hotspots different amounts of color charge by
 * applying individual weight factors. See `sample_hotspot_weights()`.
 */
class HotspotNucleus : public Nucleus {
 private:
  uint m_num_hotspots_per_nucleon;
  uint m_num_hotspots_total;
  // 1d size parameter of hotspots (std dev of normal distribution)
  double m_hotspot_size = std::sqrt(0.7);
  // array containing HotspotPos objects for all the hotspots
  HotspotPos* m_hotspot_pos = nullptr;
  // array containing hotspots weight factors for all hotspots
  double* m_hotspot_weights = nullptr;

 public:
  // set size parameter of nucleon in GeV⁻¹
  void set_nucleon_size(double nucleon_size) override;
  // set size parameter of hotspot in GeV⁻¹
  void set_hotspot_size(double hotspot_size);

  // sample only hotspot positions but leave nucleon positions untouched
  void sample_only_hotspot_pos();
  // sample the whole HotspotNucleus again (both nucleon and hotspot positions)
  void sample() override;  // also resamples the hotspot positions

  /* get nucleus thickness but using the hotspots instead of full nucleons
   * (equivalent to hotspot density)
   * NOTE: of course, this means that integration of this property over the
   * whole transverse plane will give the number of hotspots
   */
  double get_hotspot_thickness(double x, double y) const;
  uint get_num_hotspots_per_nucleon() const;
  uint get_num_hotspots_total() const;

  /* return the transvese position of hotspot `hotspot_num` in nucleon
   * `nucleon_num` -both paramters are zero-indexed, of course:
   * --`hotspot_num` is relative to a single nucleon → `hotspot_num = 0, ...,
   * num_hotspots_per_nucleon-1`
   * --`nucleon_num = 0, ..., atomic_num-1`
   */
  const HotspotPos* get_hotspot_pos(uint nucleon_num, uint hotspot_num) const;
  /* get the transverse position of the nucleon with index
   * `hotspot_num_absolute`
   * -this is zero-indexed as well!
   * -example: assuming we have 5 nucleons with 3 hotspots each:
   * --then `hotspot_num_absolute = 9` will refer to the 0th hotspot in nucleon
   * 3 (again, both are zero-indexed)
   * -this function is slightly faster than the one with two parameters because
   * we can access the array without multiplication
   */
  const HotspotPos* get_hotspot_pos(uint hotspot_num_absolute) const;

  /* sample the hotspot weights (optional)
   * -these will be multiplied onto any previously sampled nucleon weights
   * --meaning: if the respective nucleon already has a weight of `2.0` and the
   * one hotspot gets a weight of `1.5`, then the actual weight of the hotspot
   * will be set to `3.0`
   */
  void sample_hotspot_weights();
  /* sample the hotspot weights (optional) but with fixed average
   * -same as `samples_hotspot_weights()` but the individual weights will be
   * rescaled such that the total average of weights is `1.0`
   */
  void sample_hotspot_weights_fixed_avg();
  // reset the hotspot weights to `1.0`
  void reset_hotspot_weights();
  // reset all weights (nucleon and hotspot), this just calls
  // `reset_nucleon_weights()` and `reset_hotspot_weights()`
  void reset_weights();
  // get hotspot weight of hotspot `hotspot_num` in nucleon `nucleon_num`
  double get_hotspot_weight(uint nucleon_num, uint hotspot_num) const;
  // get hotspot weight of hotspot with absolute index `hotspot_num_absoulute`
  // (see `get_hotspot_pos(uint)` and `get_hotspot_pos(uint, uint)` for
  // explanations of how the paramters work in detail)
  double get_hotspot_weight(uint hotspot_num_absolute) const;

  // delete default constructor because we need to know the number of hotspots
  // per nucleon in order to allocate the correct amount of memory
  HotspotNucleus() = delete;
  // constructor will sample the base Nucleus class, then sample hotspot
  // positions for the HotspotNucleus for the specified numbers of nucleons and
  // hotspots
  HotspotNucleus(uint seed, uint atomic_num = 1,
                 uint num_hotspots_per_nucleon = 3,
                 double nucleon_size = std::sqrt(3.3),
                 double hotspot_size = std::sqrt(0.7),
                 SamplingDistribution sampling_distribution = WoodsSaxon);
  // copy constructor
  HotspotNucleus(const HotspotNucleus&);
  // move constructor
  HotspotNucleus(HotspotNucleus&&);
  // copy assignment
  HotspotNucleus& operator=(const HotspotNucleus&);
  // move assignment
  HotspotNucleus& operator=(HotspotNucleus&&);
  // destructor
  ~HotspotNucleus();

 private:
  // normal distribution for sampling of hotspot positions (nucleon
  // with 2d-Gaussian profile)
  std::normal_distribution<double> m_dist_gaussian =
      std::normal_distribution<double>(0.0, m_nucleon_size);

  // the random function returning numbers from the normal dist
  inline double m_rand_gaussian() { return m_dist_gaussian(*m_rng); }

  void safe_delete_hotspot_pos();
  void prepare_hotspot_pos();
  void safe_delete_hotspot_weights();
  void prepare_hotspot_weights();
  void sample_hotspots_single_nucleon(uint nucleon_num);
  void sample_single_hotspot_pos(HotspotPos* hotspot_pos);
  bool fits_hotspot_distribution(double r_sqr);
};

#endif  // NUCLEUS_HOTSPOTNUCLEUS_HPP_
