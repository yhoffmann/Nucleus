#pragma once


#include <stdlib.h>
#include <string>
#include <math.h>
#include <random>
#include "constants.hpp"


enum class SamplingDistribution : unsigned char
{
    WoodsSaxon,
    Gaussian
    //ShellModelWaveFunctions
};


class Nucleus
{
protected:

    uint m_atomic_num; // atomic number
    double m_mean_bulk_radius; // GeVm1 // avg radius of nuclei
    double m_mean_surface_diffusiveness; // GeVm1 // nucleus surface diffusiveness
    double* m_nucleon_pos = nullptr; // 3D positions of nucleons, relative to center of mass
    double m_nucleon_size = 3.3; // GeVm2 // interaction cross section of individual nucleon // TODO

    double m_sampling_range;
    SamplingDistribution m_sampling_distribution;
    std::mt19937& m_rng;
    std::uniform_real_distribution<double> m_dist_01u = std::uniform_real_distribution<double>(0.0, 1.0);
    inline double m_rand() { return m_dist_01u(m_rng); }

    void set_mean_bulk_radius();
    void set_mean_nucleus_diffusiveness();
    void set_sampling_range();
    void safe_delete_pos();
    void prepare_pos();
    void sample_single_pos(double pos[3]);
    virtual bool fits_distribution(double r_sqr);

public:

    void sample_nucleon_pos();
    void export_nucleon_positions(const double impact_param[2], const std::string& filepath) const;
    virtual double get_nucleus_thickness(double x, double y) const;
    uint get_atomic_num() const;
    const double* get_nucleon_pos(uint nucleon_num) const;
    void safe_get_nucleon_pos(double pos[3], uint nucleon_num) const;
    double get_mean_bulk_radius() const;
    double get_mean_surface_diffusiveness() const;
    void set_nucleon_size(double sigma_nn);
    double get_nucleon_size() const;

    Nucleus() = delete;
    Nucleus(uint atomic_num, std::mt19937& rng, SamplingDistribution sampling_distribution = SamplingDistribution::WoodsSaxon);
    Nucleus(const Nucleus& other);
    Nucleus(Nucleus&&) = delete;
    Nucleus& operator=(const Nucleus& other);
    Nucleus& operator=(Nucleus&& other);
    ~Nucleus();
};


namespace NormedSamplingDistributions
{
    inline double woods_saxon (double r, double mean_bulk_radius, double mean_surface_diffusiveness)
    {
        return ( 1.0+exp(mean_bulk_radius/mean_surface_diffusiveness) ) / ( 1.0+exp( (r-mean_bulk_radius)/mean_surface_diffusiveness ) );
    }


    inline double gaussian (double r_sqr, double sigma_sqr)
    {
        return exp( -r_sqr/(2.0*sigma_sqr) );
    }
}