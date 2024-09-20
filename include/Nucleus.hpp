#pragma once


#include <stdlib.h>
#include <string>
#include <math.h>
#include <random>


enum SamplingDistribution : unsigned char
{
    WoodsSaxon,
    Gaussian
    //ShellModelWaveFunctions
};


struct NucleonPos
{
    double x, y, z;

    friend std::ostream& operator<<(std::ostream& stream, const NucleonPos& pos);
};


class Nucleus
{
protected:

    uint m_atomic_num; // atomic number
    double m_mean_bulk_radius; // GeVm1 // avg radius of nuclei
    double m_mean_surface_diffusiveness; // GeVm1 // nucleus surface diffusiveness
    NucleonPos* m_nucleon_pos = nullptr; // 3D positions of nucleons, relative to center of mass
    double m_nucleon_size = std::sqrt(3.3); // GeVm1 // 1d size parameter of nucleon (for example std dev of normal distribution)

    double m_sampling_range;
    SamplingDistribution m_sampling_distribution;
    std::mt19937* m_rng = nullptr;
    std::uniform_real_distribution<double> m_dist_01u = std::uniform_real_distribution<double>(0.0, 1.0);
    
    inline double m_rand() { return m_dist_01u(*m_rng); }

public:

    virtual void sample_nucleon_pos();
    void export_nucleon_positions(double impact_param_x, double impact_param_y, const std::string& filepath) const;
    double get_nucleon_thickness(double x, double y) const;
    uint get_atomic_num() const;
    const NucleonPos* get_nucleon_pos(uint nucleon_num) const;
    const NucleonPos* safe_get_nucleon_pos(uint nucleon_num) const;
    double get_mean_bulk_radius() const;
    double get_mean_surface_diffusiveness() const;
    virtual void set_nucleon_size(double sigma_nn);
    double get_nucleon_size() const;

    Nucleus() = delete;
    Nucleus(uint atomic_num, uint seed, SamplingDistribution sampling_distribution = WoodsSaxon);
    Nucleus(const Nucleus&);
    Nucleus(Nucleus&&);
    Nucleus& operator=(const Nucleus&);
    Nucleus& operator=(Nucleus&&);
    virtual ~Nucleus();

protected:

    void set_mean_bulk_radius();
    void set_mean_surface_diffusiveness();
    void set_sampling_range();
    void safe_delete_pos();
    void prepare_pos();
    void safe_delete_rng();
    void prepare_rng(uint seed);
    void prepare_rng(const std::mt19937& rng);
    void sample_single_pos(NucleonPos* nucleon_pos);
    bool fits_nucleon_distribution(double r_sqr);
};


namespace NormedSamplingDistributions
{
    inline double woods_saxon (double r, double mean_bulk_radius, double mean_surface_diffusiveness)
    {
        return ( 1.0+exp(-mean_bulk_radius/mean_surface_diffusiveness) ) / ( 1.0+exp( (r-mean_bulk_radius)/mean_surface_diffusiveness ) );
    }


    inline double gaussian (double r_sqr, double sigma_sqr)
    {
        return exp( -r_sqr/(2.0*sigma_sqr) );
    }
}