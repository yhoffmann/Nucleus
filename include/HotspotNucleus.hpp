#pragma once
#ifndef HOTSPOT_NUCLEUS_HPP
#define HOTSPOT_NUCLEUS_HPP


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
    uint m_num_hotspots_total = 3;
    double m_hotspot_size = std::sqrt(0.7); // 1d size parameter of hotspots (std dev of normal distribution)
    HotspotPos* m_hotspot_pos = nullptr;
    double* m_hotspot_weights = nullptr;

public:

    void set_nucleon_size(double nucleon_size) override;
    void set_hotspot_size(double hotspot_size);

    void sample_only_hotspot_pos();
    void sample() override; // also resamples the hotspot positions

    double get_hotspot_thickness(double x, double y) const;
    uint get_num_hotspots_per_nucleon() const;
    uint get_num_hotspots_total() const;

    const HotspotPos* get_hotspot_pos(uint nucleon_num, uint hotspot_num) const;
    const HotspotPos* get_hotspot_pos(uint hotspot_num_absolute) const; // provide absolute hotspot number; e.g. 4 if you want to access hotspot 1 in nucleon 1 (both 0-indexed!!!) in the case of 3 hotspots per nucleon and at least 2 nucleons

    void sample_hotspot_weights();
    void sample_hotspot_weights_fixed_avg();
    double get_hotspot_weight(uint nucleon_num, uint hotspot_num) const;
    double get_hotspot_weight(uint hotspot_num_absolute) const;

    HotspotNucleus() = delete;
    HotspotNucleus(uint seed, uint atomic_num = 1, uint num_hotspots_per_nucleon = 3, double nucleon_size = std::sqrt(3.3), double hotspot_size = std::sqrt(0.7), SamplingDistribution sampling_distribution = WoodsSaxon);
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
    void safe_delete_hotspot_weights();
    void prepare_hotspot_weights();
    void sample_hotspots_single_nucleon(uint nucleon_num);
    void sample_single_hotspot_pos(HotspotPos* hotspot_pos);
    bool fits_hotspot_distribution(double r_sqr);
};


#endif // HOTSPOT_NUCLEUS_HPP