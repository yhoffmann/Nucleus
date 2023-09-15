#pragma once


#include <stdlib.h>
#include <string>
#include <math.h>


class NNEvent;


enum class SamplingDistribution : unsigned char
{
    WoodsSaxon,
    ShellModelWaveFunctions
};


class Nucleus
{
    uint m_atomic_num; // atomic number
    double m_bulk_radius; // avg radius of nuclei
    double m_mean_surface_diffusiveness; // nucleus surface diffusiveness
    double* m_pos = nullptr; // 3D positions of nucleons, relative to center of mass
    double m_sigma_nn; // interaction cross section of individual nucleon
    double m_sqrt_s_nn; // com energy of collision (per nucleon) // TODO check if this is true

    void set_bulk_radius(uint atomic_num);
    void set_mean_nucleus_diffusiveness(uint atomic_num);
    void set_sigma_nn();
    void set_sqrt_s_nn();
    void safe_delete_pos();
    void prepare_pos();
    void sample_single_pos(double pos[3]);
    bool check_fits_distribution(double param, SamplingDistribution dist) const;
    void get_nucleon_pos(double pos[3], uint nucleon_num) const;

public:

    void sample_pos();
    void export_nucleon_positions(const double impact_param[2], const std::string& filepath) const;
    double get_nucleus_thickness(double x, double y) const;
    uint get_atomic_num() const;
    void safe_get_nucleon_pos(double pos[3], uint nucleon_num) const;
    double get_bulk_radius() const;
    double get_mean_surface_diffusiveness() const;
    double get_sigma_nn() const;
    double get_sqrt_s_nn() const;

    Nucleus() = delete;
    Nucleus(uint atomic_num);
    Nucleus(const Nucleus& other);
    Nucleus(Nucleus&&) = delete;
    Nucleus& operator=(const Nucleus& other);
    Nucleus& operator=(Nucleus&& other);
    ~Nucleus();

    friend class NNEvent;
};


namespace SamplingDistributions
{
    double woods_saxon_no_rho_0_inverse_cdf(double rng_param, double bulk_radius, double mean_surface_diffusiveness);
}