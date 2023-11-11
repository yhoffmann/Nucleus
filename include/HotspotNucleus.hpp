#pragma once


#include <random>
#include "Nucleus.hpp"


class HotspotNucleus : public Nucleus
{
public:

    void set_nucleon_size(double nucleon_size) override;
    void set_hotspot_size(double hotspot_size);
    void sample_hotspot_pos();
    void sample_nucleon_pos() override; // also resamples the hotspot positions
    double get_hotspot_thickness(double x, double y) const;
    uint get_num_hotspots_per_nucleon() const;
    const double* get_hotspot_pos(uint nucleon_num, uint hotspot_num) const;

    HotspotNucleus() = delete;
    HotspotNucleus(uint atomic_num, uint num_hotspots_per_nucleon, std::mt19937& rng, SamplingDistribution sampling_distribution = SamplingDistribution::WoodsSaxon);
    HotspotNucleus(const HotspotNucleus&) = delete;
    HotspotNucleus(HotspotNucleus&&) = delete;
    ~HotspotNucleus();

private:

    uint m_num_hotspots_per_nucleon;
    double* m_hotspot_pos = nullptr;
    double m_hotspot_size = std::sqrt(0.7); // 1d size parameter of hotspots (std dev of normal distribution)

    std::normal_distribution<double> m_dist_gaussian = std::normal_distribution<double>(0.0, m_nucleon_size);
    inline double m_rand_gaussian() { return m_dist_gaussian(m_rng); }

    void safe_delete_hotspot_pos();
    void prepare_hotspot_pos();
    void sample_hotspots_single_nucleon(uint nucleon_num);
    void sample_single_hotspot_pos(double hotspot_pos[2]);
    bool fits_hotspot_distribution(double r_sqr);
};