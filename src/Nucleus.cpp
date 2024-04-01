#include <iostream>
#include <math.h>
#include <fstream>
#include <cstring>
#include <random>
#include <stdlib.h>
#include "../include/Nucleus.hpp"
#include "../include/constants.hpp"
#include "../include/NuclearParameters.hpp"


std::ostream& operator<<(std::ostream& stream, const NucleonPos& pos)
{
    stream << pos.x<<" "<<pos.y<<" "<<pos.z;

    return stream;
}


void Nucleus::sample_nucleon_pos()
{
    NucleonPos center_of_mass = {0.0, 0.0, 0.0};

    for (uint n=0; n<m_atomic_num; n++)
    {
        sample_single_pos(m_nucleon_pos+n);

        center_of_mass.x += m_nucleon_pos[n].x;
        center_of_mass.y += m_nucleon_pos[n].y;
        center_of_mass.z += m_nucleon_pos[n].z;
    }

    double inverse_divisor = 1.0/double(m_atomic_num);
    center_of_mass.x *= inverse_divisor;
    center_of_mass.y *= inverse_divisor;
    center_of_mass.z *= inverse_divisor;

    for (uint n=0; n<m_atomic_num; n++)
    {
        m_nucleon_pos[n].x -= center_of_mass.x;
        m_nucleon_pos[n].y -= center_of_mass.y;
        m_nucleon_pos[n].z -= center_of_mass.z;
    }
}


void Nucleus::export_nucleon_positions (double impact_param_x, double impact_param_y, const std::string& filename) const
{
    std::ofstream filestream;
    filestream.open(filename);
    if (!filestream.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open file " << filename << std::endl;
#endif
        exit(33);
    }

    filestream << "# " << m_mean_bulk_radius << " " << impact_param_x << " " << impact_param_y << std::endl; 

    for (uint n=0; n<m_atomic_num; n++)
    {
        filestream << m_nucleon_pos[n].x << " " << m_nucleon_pos[n].y << " " << m_nucleon_pos[n].z << " " << impact_param_x << " " << impact_param_y << " " << m_nucleon_size << std::endl;
    }

    filestream.close();
}


double Nucleus::get_nucleon_thickness (double x, double y) const
{
    double thickness = 0.0;

    double inverse_divisor = 1.0/(2.0*m_nucleon_size*m_nucleon_size);
    for (uint n=0; n<m_atomic_num; n++)
    {
        double delta_x = x-m_nucleon_pos[n].x;
        double delta_y = y-m_nucleon_pos[n].y;

        double r_sqr = delta_x*delta_x + delta_y*delta_y;

        thickness += exp( -r_sqr*inverse_divisor );
    }
    return thickness*inverse_divisor/M_PI;
}


uint Nucleus::get_atomic_num() const
{
    return m_atomic_num;
}


const NucleonPos* Nucleus::get_nucleon_pos (uint nucleon_num) const
{
    return m_nucleon_pos+nucleon_num;
}


const NucleonPos* Nucleus::safe_get_nucleon_pos (uint nucleon_num) const
{
    if (nucleon_num>m_atomic_num-1)
    {
#ifndef _QUIET
        std::cerr << "Accessing nucleon outside of range. Truncating to last nucleon..." << std::endl;
#endif
        nucleon_num = m_atomic_num-1;
    }

    return get_nucleon_pos(nucleon_num);
}


double Nucleus::get_mean_bulk_radius() const
{
    return m_mean_bulk_radius;
}


double Nucleus::get_mean_surface_diffusiveness() const
{
    return m_mean_surface_diffusiveness;
}


double Nucleus::get_nucleon_size() const
{
    return m_nucleon_size;
}


void Nucleus::set_nucleon_size (double nucleon_size)
{
    m_nucleon_size = nucleon_size;
}


Nucleus::Nucleus (uint atomic_num, std::mt19937& rng, SamplingDistribution sampling_distribution)
    : m_atomic_num(atomic_num)
    , m_sampling_distribution(sampling_distribution)
    , m_rng(rng)
{
    set_mean_bulk_radius();
    set_mean_surface_diffusiveness();
    
    set_sampling_range();

    prepare_pos();
    sample_nucleon_pos();
}


Nucleus::Nucleus (const Nucleus& other)
    : m_atomic_num(other.m_atomic_num)
    , m_mean_bulk_radius(other.m_mean_bulk_radius)
    , m_mean_surface_diffusiveness(other.m_mean_surface_diffusiveness)
    , m_nucleon_size(other.m_nucleon_size)
    , m_sampling_range(other.m_sampling_range)
    , m_sampling_distribution(other.m_sampling_distribution)
    , m_rng(other.m_rng)
{
    prepare_pos();
    std::copy(other.m_nucleon_pos, other.m_nucleon_pos+3*m_atomic_num, m_nucleon_pos);
}


Nucleus::Nucleus (Nucleus&& other)
    : m_atomic_num(other.m_atomic_num)
    , m_mean_bulk_radius(other.m_mean_bulk_radius)
    , m_mean_surface_diffusiveness(other.m_mean_surface_diffusiveness)
    , m_nucleon_size(other.m_nucleon_size)
    , m_sampling_range(other.m_sampling_range)
    , m_sampling_distribution(other.m_sampling_distribution)
    , m_rng(other.m_rng)
{
    safe_delete_pos();
    m_nucleon_pos = other.m_nucleon_pos;

    other.m_nucleon_pos = nullptr;
}


Nucleus& Nucleus::operator= (const Nucleus& other)
{
    if (this==&other)
        return *this;

    safe_delete_pos();

    m_atomic_num = other.m_atomic_num;
    m_mean_bulk_radius = other.m_mean_bulk_radius;
    m_mean_surface_diffusiveness = other.m_mean_surface_diffusiveness;
    m_nucleon_size = other.m_nucleon_size;
    
    m_sampling_range = other.m_sampling_range;
    m_sampling_distribution = other.m_sampling_distribution;
    m_rng = other.m_rng;

    prepare_pos();
    std::copy(other.m_nucleon_pos, other.m_nucleon_pos+m_atomic_num, m_nucleon_pos);

    return *this;
}


Nucleus& Nucleus::operator= (Nucleus&& other)
{
    if (this==&other)
        return *this;

    safe_delete_pos();

    m_atomic_num = other.m_atomic_num;
    m_mean_bulk_radius = other.m_mean_bulk_radius;
    m_mean_surface_diffusiveness = other.m_mean_surface_diffusiveness;
    m_nucleon_size = other.m_nucleon_size;
    
    m_sampling_range = other.m_sampling_range;
    m_sampling_distribution = other.m_sampling_distribution;
    m_rng = other.m_rng;

    other.m_nucleon_pos = nullptr;

    return *this;
}


Nucleus::~Nucleus()
{
    safe_delete_pos();
}


void Nucleus::set_mean_bulk_radius()
{
    m_mean_bulk_radius = fmToGeVm1*NuclearParameters::mean_bulk_radius(m_atomic_num);
}


void Nucleus::set_mean_surface_diffusiveness()
{
    m_mean_surface_diffusiveness = fmToGeVm1*NuclearParameters::mean_surface_diffusiveness(m_atomic_num);
}


void Nucleus::set_sampling_range()
{
    switch (m_sampling_distribution)
    {
        case SamplingDistribution::WoodsSaxon:
            m_sampling_range = m_mean_bulk_radius+10.0*m_mean_surface_diffusiveness;
        break;

        case SamplingDistribution::Gaussian:
        break;
    }
}


void Nucleus::safe_delete_pos()
{
    if (m_nucleon_pos)
        delete[] m_nucleon_pos;
    m_nucleon_pos = nullptr;
}


void Nucleus::prepare_pos()
{
    safe_delete_pos();

    m_nucleon_pos = new(std::nothrow) NucleonPos [m_atomic_num];
    if (m_nucleon_pos == nullptr)
        exit(31);
}


void Nucleus::sample_single_pos (NucleonPos* nucleon_pos)
{
    while(true)
    {
        NucleonPos sampled_pos;

        sampled_pos.x = m_sampling_range*2.0*(m_rand()-0.5);
        sampled_pos.y = m_sampling_range*2.0*(m_rand()-0.5);
        sampled_pos.z = m_sampling_range*2.0*(m_rand()-0.5);

        double r_sqr = sampled_pos.x*sampled_pos.x + sampled_pos.y*sampled_pos.y + sampled_pos.z*sampled_pos.z;

        if (fits_nucleon_distribution(r_sqr))
        {
            nucleon_pos->x = sampled_pos.x;
            nucleon_pos->y = sampled_pos.y;
            nucleon_pos->z = sampled_pos.z;

            return;
        }
    }
}


bool Nucleus::fits_nucleon_distribution (double r_sqr)
{
    switch (m_sampling_distribution)
    {
        case SamplingDistribution::WoodsSaxon:
            return m_rand() < NormedSamplingDistributions::woods_saxon(std::sqrt(r_sqr), m_mean_bulk_radius, m_mean_surface_diffusiveness);

        case SamplingDistribution::Gaussian:
            return true;

        default:
            return true;
    }
}