#include "../include/Nucleus.hpp"
#include <iostream>
#include "../include/constants.hpp"
#include <math.h>
#include <fstream>
#include <cstring>
#include <random>


class NNEvent;


void Nucleus::set_bulk_radius (uint atomic_num)
{
    if (atomic_num==208)
        m_bulk_radius = 6.62;
    else if (atomic_num==197)
        m_bulk_radius = 6.38;
    else if (atomic_num==2000)
        m_bulk_radius = 5;
    else
    {
        std::cerr << "No nuclear data for this atomic number" << std::endl;
        exit(0);
    }
}


void Nucleus::set_mean_nucleus_diffusiveness (uint atomic_num)
{
    if (atomic_num==208)
        m_mean_surface_diffusiveness = 0.546;
    else if (atomic_num==197)
        m_mean_surface_diffusiveness = 0.535;
    else if (atomic_num==2000)
        m_mean_surface_diffusiveness = 0.005;
    else
    {
        std::cerr << "No nuclear data for this atomic number" << std::endl;
        exit(0);
    }
}


void Nucleus::set_sigma_nn()
{
    if (m_atomic_num==197)
        m_sigma_nn = 42e6*nbTofm2;
    if (m_atomic_num==208)
        m_sigma_nn = 64e6*nbTofm2;
}


void Nucleus::set_sqrt_s_nn()
{
    if (m_atomic_num==197)
        m_sqrt_s_nn = 200*GeVTofmm1;
    if (m_atomic_num==208)
        m_sqrt_s_nn = 2.76e3*GeVTofmm1;
}

void Nucleus::safe_delete_pos()
{
    if (m_pos)
        delete[] m_pos;
    m_pos = nullptr;
}


void Nucleus::prepare_pos()
{
    safe_delete_pos();

    m_pos = new double [m_atomic_num*3];
}


void Nucleus::sample_single_pos (double pos[3])
{
    double r_max = 4.0*m_bulk_radius;

    bool fits_dstribution = 0;
    while(!fits_dstribution)
    {
        pos[0] = r_max*2.0*(drand48()-0.5);
        pos[1] = r_max*2.0*(drand48()-0.5);
        pos[2] = r_max*2.0*(drand48()-0.5);

        double r_sqr = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
        fits_dstribution = check_fits_distribution(r_sqr, SamplingDistribution::WoodsSaxon);
    }
}


bool Nucleus::check_fits_distribution(double r_sqr, SamplingDistribution dist) const
{
    if (dist==SamplingDistribution::WoodsSaxon)
        return r_sqr<std::pow(SamplingDistributions::woods_saxon_no_rho_0_inverse_cdf(drand48(), m_bulk_radius, m_mean_surface_diffusiveness), 2);
    else return true;
}


void Nucleus::sample_pos()
{
    double center_of_mass[3] = {0.0, 0.0, 0.0};

    for (uint n=0; n<m_atomic_num; n++)
    {
        sample_single_pos(m_pos+n*3);

        center_of_mass[0] += m_pos[n*3];
        center_of_mass[1] += m_pos[n*3+1];
        center_of_mass[2] += m_pos[n*3+2];
    }

    center_of_mass[0] /= double(m_atomic_num);
    center_of_mass[1] /= double(m_atomic_num);
    center_of_mass[2] /= double(m_atomic_num);

    for (uint n=0; n<m_atomic_num; n++)
    {
        m_pos[n*3] -= center_of_mass[0];
        m_pos[n*3+1] -= center_of_mass[1];
        m_pos[n*3+2] -= center_of_mass[2];
    }
}


void Nucleus::export_nucleon_positions (const double impact_parameter[2], const std::string& filename) const
{
    std::ofstream filestream;
    filestream.open("Data/"+filename);
    if (!filestream.is_open())
    {
        std::cerr << "Could not open file " << "Data/"+filename << std::endl;
        exit(0);
    }

    for (uint n=0; n<m_atomic_num; n++)
    {
        filestream << m_pos[n*3] << " " << m_pos[n*3+1] << " " << m_pos[n*3+2] << " " << impact_parameter[0] << " " << impact_parameter[1] << " " << std::sqrt(m_sigma_nn/M_PI)/2.0 << std::endl;
    }
    filestream.close();

    filestream.open("Data/Radius"+filename);
    if (!filestream.is_open())
    {
        std::cerr << "Could not open file " << "Data/Radius"+filename << std::endl;
        exit(0);
    }

    filestream << m_bulk_radius << " " << impact_parameter[0] << " " << impact_parameter[1] << std::endl;

    filestream.close();
}


double Nucleus::get_nucleus_thickness (double x, double y) const
{
    double thickness = 0.0;

    for (uint n=0; n<m_atomic_num; n++)
    {
        double r_sqr = (m_pos[n*3]-x)*(m_pos[n*3]-x)+(m_pos[n*3+1]-y)*(m_pos[n*3+1]-y);

        thickness += exp(-r_sqr/(2.0*m_sigma_nn));
    }
    return thickness/(2.0*M_PI*m_sigma_nn);
}


uint Nucleus::get_atomic_num() const
{
    return m_atomic_num;
}


const double* Nucleus::get_nucleon_pos (uint nucleon_num) const
{
    return m_pos+3*nucleon_num;
}


void Nucleus::safe_get_nucleon_pos (double pos[3], uint nucleon_num) const
{
    if (nucleon_num>m_atomic_num-1)
    {
        std::cerr << "Accessing nucleon outside of range. Truncating to last nucleon..." << std::endl;
        nucleon_num = m_atomic_num-1;
    }

    const double* pos_ptr = get_nucleon_pos(nucleon_num);

    pos[0] = *pos_ptr;
    pos[1] = *(pos_ptr+1);
    pos[2] = *(pos_ptr+2);
}


double Nucleus::get_bulk_radius() const
{
    return m_bulk_radius;
}


double Nucleus::get_mean_surface_diffusiveness() const
{
    return m_mean_surface_diffusiveness;
}


double Nucleus::get_sigma_nn() const
{
    return m_sigma_nn;
}


double Nucleus::get_sqrt_s_nn() const
{
    return m_sqrt_s_nn;
}


Nucleus::Nucleus (uint atomic_num)
    : m_atomic_num(atomic_num)
{
    set_bulk_radius(m_atomic_num);
    set_mean_nucleus_diffusiveness(m_atomic_num);

    set_sigma_nn();
    set_sqrt_s_nn();
    
    prepare_pos();

    sample_pos();
}


Nucleus::Nucleus (const Nucleus& other)
    : m_atomic_num(other.m_atomic_num), m_bulk_radius(other.m_bulk_radius), m_mean_surface_diffusiveness(other.m_mean_surface_diffusiveness), m_sigma_nn(other.m_sigma_nn), m_sqrt_s_nn(other.m_sqrt_s_nn)
{
    prepare_pos();

    memcpy(m_pos, other.m_pos, 3*m_atomic_num*sizeof(double));
}


Nucleus& Nucleus::operator= (const Nucleus& other)
{
    if (this==&other)
        return *this;

    safe_delete_pos();

    memcpy(this, &other, sizeof(Nucleus));

    m_pos = new double [3*m_atomic_num];
    memcpy(m_pos, other.m_pos, 3*m_atomic_num*sizeof(double));

    return *this;
}


Nucleus& Nucleus::operator= (Nucleus&& other)
{
    if (this==&other)
        return *this;

    safe_delete_pos();

    memcpy(this, &other, sizeof(Nucleus));

    other.m_pos = nullptr;

    return *this;
}


Nucleus::~Nucleus()
{
    safe_delete_pos();
}


namespace SamplingDistributions
{
    double woods_saxon_no_rho_0_inverse_cdf(double rng_param, double bulk_radius, double mean_surface_diffusiveness)
    {
        return mean_surface_diffusiveness*log(1.0/rng_param-1.0)+bulk_radius;
    }
}