#pragma once


#include <stdlib.h>
#include <string>
#include <math.h>
#include <random>


enum class SamplingDistribution : unsigned char
{
    WoodsSaxon,
    //ShellModelWaveFunctions
};


class Nucleus
{
    std::mt19937& m_rng;
    std::uniform_real_distribution<double> m_dist = std::uniform_real_distribution<double>(0.0, 1.0);
    inline double m_rand() { return m_dist(m_rng); }

    uint m_atomic_num; // atomic number
    double m_mean_bulk_radius; // avg radius of nuclei
    double m_mean_surface_diffusiveness; // nucleus surface diffusiveness
    double* m_pos = nullptr; // 3D positions of nucleons, relative to center of mass
    double m_sigma_nn = 0.5*0.5*M_PI; //fm2 // interaction cross section of individual nucleon // TODO

    void set_mean_bulk_radius();
    void set_mean_nucleus_diffusiveness();
    void safe_delete_pos();
    void prepare_pos();
    void sample_single_pos(double pos[3]);
    bool check_fits_distribution(double param, SamplingDistribution dist);

public:

    void sample();
    void export_nucleon_positions(const double impact_param[2], const std::string& filepath) const;
    double get_nucleus_thickness(double x, double y) const;
    uint get_atomic_num() const;
    const double* get_nucleon_pos(uint nucleon_num) const;
    void safe_get_nucleon_pos(double pos[3], uint nucleon_num) const;
    double get_mean_bulk_radius() const;
    double get_mean_surface_diffusiveness() const;
    double get_sigma_nn() const;
    double get_sqrt_s_nn() const;

    Nucleus() = delete;
    Nucleus(std::mt19937& rng, uint atomic_num);
    Nucleus(const Nucleus& other);
    Nucleus(Nucleus&&) = delete;
    Nucleus& operator=(const Nucleus& other);
    Nucleus& operator=(Nucleus&& other);
    ~Nucleus();
};


namespace SamplingDistributions
{
    double woods_saxon_no_rho_0_inverse_cdf(double rng_param, double bulk_radius, double mean_surface_diffusiveness);
}