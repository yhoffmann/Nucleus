#include <iostream>
#include <math.h>
#include <fstream>
#include <cstring>
#include <random>
#include <stdlib.h>
#include "../include/Nucleus.hpp"
#include "../include/constants.hpp"
#include "../include/NuclearParameters.hpp"


void Nucleus::sample_nucleon_pos()
{
    double center_of_mass[3] = {0.0, 0.0, 0.0};

    for (uint n=0; n<m_atomic_num; n++)
    {
        sample_single_pos(m_nucleon_pos+n*3);

        center_of_mass[0] += m_nucleon_pos[n*3];
        center_of_mass[1] += m_nucleon_pos[n*3+1];
        center_of_mass[2] += m_nucleon_pos[n*3+2];
    }

    center_of_mass[0] /= double(m_atomic_num);
    center_of_mass[1] /= double(m_atomic_num);
    center_of_mass[2] /= double(m_atomic_num);

    for (uint n=0; n<m_atomic_num; n++)
    {
        m_nucleon_pos[n*3] -= center_of_mass[0];
        m_nucleon_pos[n*3+1] -= center_of_mass[1];
        m_nucleon_pos[n*3+2] -= center_of_mass[2];
    }
}


void Nucleus::export_nucleon_positions (const double impact_parameter[2], const std::string& filename) const
{
    std::ofstream filestream;
    filestream.open("Data/"+filename);
    if (!filestream.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open file " << "Data/"+filename << std::endl;
#endif
        exit(33);
    }

    for (uint n=0; n<m_atomic_num; n++)
    {
        filestream << m_nucleon_pos[n*3] << " " << m_nucleon_pos[n*3+1] << " " << m_nucleon_pos[n*3+2] << " " << impact_parameter[0] << " " << impact_parameter[1] << " " << std::sqrt(m_nucleon_size/M_PI)/2.0 << std::endl; // TODO fix size
    }
    filestream.close();

    filestream.open("Data/Radius"+filename);
    if (!filestream.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open file " << "Data/Radius"+filename << std::endl;
#endif
        exit(33);
    }

    filestream << m_mean_bulk_radius << " " << impact_parameter[0] << " " << impact_parameter[1] << std::endl;

    filestream.close();
}


double Nucleus::get_nucleon_thickness (double x, double y) const
{
    double thickness = 0.0;

    for (uint n=0; n<m_atomic_num; n++)
    {
        double delta_x = x-m_nucleon_pos[n*3];
        double delta_y = y-m_nucleon_pos[n*3+1];

        double r_sqr = delta_x*delta_x + delta_y*delta_y;

        thickness += exp( -r_sqr/(2.0*m_nucleon_size) );
    }
    return thickness/(2.0*M_PI*m_nucleon_size);
}


uint Nucleus::get_atomic_num() const
{
    return m_atomic_num;
}


const double* Nucleus::get_nucleon_pos (uint nucleon_num) const
{
    return m_nucleon_pos+3*nucleon_num;
}


void Nucleus::safe_get_nucleon_pos (double pos[3], uint nucleon_num) const
{
    if (nucleon_num>m_atomic_num-1)
    {
#ifndef _QUIET
        std::cerr << "Accessing nucleon outside of range. Truncating to last nucleon..." << std::endl;
#endif
        nucleon_num = m_atomic_num-1;
    }

    const double* pos_ptr = get_nucleon_pos(nucleon_num);

    pos[0] = pos_ptr[0];
    pos[1] = pos_ptr[1];
    pos[2] = pos_ptr[2];
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
    : m_atomic_num(atomic_num), m_sampling_distribution(sampling_distribution), m_rng(rng)
{
    set_mean_bulk_radius();
    set_mean_surface_diffusiveness();
    
    set_sampling_range();

    prepare_pos();
    sample_nucleon_pos();
}


Nucleus::Nucleus (const Nucleus& other)
    : m_atomic_num(other.m_atomic_num), m_mean_bulk_radius(other.m_mean_bulk_radius), m_mean_surface_diffusiveness(other.m_mean_surface_diffusiveness), m_nucleon_size(other.m_nucleon_size), m_sampling_range(other.m_sampling_range), m_sampling_distribution(other.m_sampling_distribution), m_rng(other.m_rng)
{
    prepare_pos();
    std::copy(other.m_nucleon_pos, other.m_nucleon_pos+3*m_atomic_num, m_nucleon_pos);
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
    std::copy(other.m_nucleon_pos, other.m_nucleon_pos+3*m_atomic_num, m_nucleon_pos);

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

    m_nucleon_pos = new double [m_atomic_num*3];
    if (m_nucleon_pos == nullptr)
        exit(31);
}


void Nucleus::sample_single_pos (double pos[3])
{
    while(true)
    {
        double sampled_pos[3];

        sampled_pos[0] = m_sampling_range*2.0*(m_rand()-0.5);
        sampled_pos[1] = m_sampling_range*2.0*(m_rand()-0.5);
        sampled_pos[2] = m_sampling_range*2.0*(m_rand()-0.5);

        double r_sqr = sampled_pos[0]*sampled_pos[0]+sampled_pos[1]*sampled_pos[1]+sampled_pos[2]*sampled_pos[2];

        if ( fits_nucleon_distribution(r_sqr) )
        {
            pos[0] = sampled_pos[0];
            pos[1] = sampled_pos[1];
            pos[2] = sampled_pos[2];

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