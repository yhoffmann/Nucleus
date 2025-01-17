#pragma once


#include <random>
#include <iostream>
#include "Nucleus.hpp"


struct HotspotPos
{
    double x, y;

    friend std::ostream& operator<<(std::ostream& stream, const HotspotPos& pos);
};


class HotspotNucleus : public Nucleus
{
private:

    uint m_num_hotspots_per_nucleon = 3;
    double m_hotspot_size = std::sqrt(0.7); // 1d size parameter of hotspots (std dev of normal distribution)
    HotspotPos* m_hotspot_pos = nullptr;

public:

    void set_nucleon_size(double nucleon_size) override;
    void set_hotspot_size(double hotspot_size);
    void sample_only_hotspot_pos();
    void sample() override; // also resamples the hotspot positions
    double get_hotspot_thickness(double x, double y) const;
    uint get_num_hotspots_per_nucleon() const;
    const HotspotPos* get_hotspot_pos(uint nucleon_num, uint hotspot_num) const;

    HotspotNucleus() = delete;
    HotspotNucleus(uint seed, uint atomic_num = 1, uint num_hotspots_per_nucleon = 3, double nucleon_size = std::sqrt(3.3), double hotspot_size = std::sqrt(0.7), double mean_bulk_radius = 0.0001, double mean_surface_diffusiveness = 0.0, SamplingDistribution sampling_distribution = WoodsSaxon);
    HotspotNucleus(const HotspotNucleus&);
    HotspotNucleus(HotspotNucleus&&);
    HotspotNucleus& operator=(const HotspotNucleus&);
    HotspotNucleus& operator=(HotspotNucleus&&);
    ~HotspotNucleus();

private:

    std::normal_distribution<double> m_dist_gaussian = std::normal_distribution<double>(0.0, m_nucleon_size);
    
    inline double m_rand_gaussian() { return m_dist_gaussian(*m_rng); }

    void safe_delete_hotspot_pos();
    void prepare_hotspot_pos();
    void sample_hotspots_single_nucleon(uint nucleon_num);
    void sample_single_hotspot_pos(HotspotPos* hotspot_pos);
    bool fits_hotspot_distribution(double r_sqr);
};