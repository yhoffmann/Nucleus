#include "../include/HotspotNucleus.hpp"


std::ostream& operator<<(std::ostream& stream, const HotspotPos& pos)
{
    stream << pos.x<<" "<<pos.y;

    return stream;
}


void HotspotNucleus::set_nucleon_size (double nucleon_size)
{
    Nucleus::set_nucleon_size(nucleon_size);
    m_dist_gaussian = std::normal_distribution<double>(0.0, nucleon_size);
}


void HotspotNucleus::set_hotspot_size (double hotspot_size)
{
    m_hotspot_size = hotspot_size;
}


void HotspotNucleus::sample_only_hotspot_pos()
{
    for (uint n=0; n<m_atomic_num; n++)
    {
        sample_hotspots_single_nucleon(n);
    }
}


void HotspotNucleus::sample()
{
    Nucleus::sample();

    sample_only_hotspot_pos();
}


double HotspotNucleus::get_hotspot_thickness (double x, double y) const
{
    double thickness = 0.0;

    double inverse_r_sqr_divisor = 1.0/(2.0*m_hotspot_size*m_hotspot_size);
    for (uint i=0; i<m_num_hotspots_total; ++i)
    {
        double delta_x = x - m_hotspot_pos[i].x;
        double delta_y = y - m_hotspot_pos[i].y;

        double r_sqr = delta_x*delta_x + delta_y*delta_y;

        thickness += exp( -r_sqr*inverse_r_sqr_divisor );
    }

    return thickness/double(m_num_hotspots_per_nucleon*M_PI)*inverse_r_sqr_divisor;
}


uint HotspotNucleus::get_num_hotspots_per_nucleon() const
{
    return m_num_hotspots_per_nucleon;
}


uint HotspotNucleus::get_num_hotspots_total() const
{
    return m_num_hotspots_total;
}


const HotspotPos* HotspotNucleus::get_hotspot_pos (uint nucleon_num, uint hotspot_num) const
{
    return &m_hotspot_pos[nucleon_num*m_num_hotspots_per_nucleon + hotspot_num];
}


const HotspotPos* HotspotNucleus::get_hotspot_pos (uint hotspot_num_absolute) const
{
    return &m_hotspot_pos[hotspot_num_absolute];
}


void HotspotNucleus::sample_hotspot_weights()
{
    prepare_hotspot_weights();

    if (nullptr == m_nucleon_weights) // nucleon weights not set
    {
        for (uint h=0; h<m_num_hotspots_total; h++)
        {
            m_hotspot_weights[h] = m_lognorm_rand();
        }
    }
    else // nucleon weights set, need to consider as well
    {
        for (uint n=0; n<m_atomic_num; n++)
        {
            for (uint h=0; h<m_num_hotspots_per_nucleon; h++)
            {
                m_hotspot_weights[n*m_num_hotspots_per_nucleon + h] = m_lognorm_rand() * get_nucleon_weight(n);
            }
        }
    }
}


double HotspotNucleus::get_hotspot_weight (uint nucleon_num, uint hotspot_num) const
{
    return m_hotspot_weights[nucleon_num*m_num_hotspots_per_nucleon + hotspot_num];
}


double HotspotNucleus::get_hotspot_weight (uint hotspot_num_absolute) const
{
    return m_hotspot_weights[hotspot_num_absolute];
}


HotspotNucleus::HotspotNucleus (uint seed, uint atomic_num, uint num_hotspots_per_nucleon, double nucleon_size, double hotspot_size, double mean_bulk_radius, double mean_surface_diffusiveness, SamplingDistribution sampling_distribution)    
    : Nucleus(seed, atomic_num, nucleon_size, mean_bulk_radius, mean_surface_diffusiveness, sampling_distribution)
    , m_num_hotspots_per_nucleon(num_hotspots_per_nucleon)
    , m_num_hotspots_total(atomic_num * num_hotspots_per_nucleon)
    , m_hotspot_size(hotspot_size)
{   
    prepare_hotspot_pos();
    sample_only_hotspot_pos();
}   


HotspotNucleus::HotspotNucleus (const HotspotNucleus& other)
    : Nucleus(other)
    , m_num_hotspots_per_nucleon(other.m_num_hotspots_per_nucleon)
    , m_num_hotspots_total(other.m_num_hotspots_total)
    , m_hotspot_size(other.m_hotspot_size)
{
    prepare_hotspot_pos();
    std::copy(other.m_hotspot_pos, other.m_hotspot_pos + m_num_hotspots_total, m_hotspot_pos);

    if (nullptr != other.m_hotspot_weights)
    {
        prepare_hotspot_weights();
        std::copy(other.m_hotspot_weights, other.m_hotspot_weights + m_num_hotspots_total, m_hotspot_weights);
    }
}


HotspotNucleus::HotspotNucleus (HotspotNucleus&& other)
    : Nucleus(std::move(other))
    , m_num_hotspots_per_nucleon(other.m_num_hotspots_per_nucleon)
    , m_num_hotspots_total(other.m_num_hotspots_total)
    , m_hotspot_size(other.m_hotspot_size)
    , m_hotspot_pos(other.m_hotspot_pos)
    , m_hotspot_weights(other.m_hotspot_weights)
{
    other.m_hotspot_pos = nullptr;
    other.m_hotspot_weights = nullptr;
}


HotspotNucleus& HotspotNucleus::operator= (const HotspotNucleus& other)
{
    if (this == &other)
        return *this;

    Nucleus::operator=(other);
    m_num_hotspots_per_nucleon = other.m_num_hotspots_per_nucleon;
    m_num_hotspots_total = other.m_num_hotspots_total;
    m_hotspot_pos = other.m_hotspot_pos;
    m_hotspot_size = other.m_hotspot_size;

    prepare_hotspot_pos();
    std::copy(other.m_hotspot_pos, other.m_hotspot_pos + m_num_hotspots_total, m_hotspot_pos);

    if (nullptr != other.m_hotspot_weights)
    {
        prepare_hotspot_weights();
        std::copy(other.m_hotspot_weights, other.m_hotspot_weights + m_num_hotspots_total, m_hotspot_weights);
    }

    return *this;
}


HotspotNucleus& HotspotNucleus::operator= (HotspotNucleus&& other)
{
    if (this == &other)
        return *this;

    Nucleus::operator=(std::move(other));
    m_num_hotspots_per_nucleon = other.m_num_hotspots_per_nucleon;
    m_num_hotspots_total = other.m_num_hotspots_total;
    m_hotspot_size = other.m_hotspot_size;
    m_hotspot_pos = other.m_hotspot_pos;
    m_hotspot_weights = other.m_hotspot_weights;

    other.m_hotspot_pos = nullptr;
    other.m_hotspot_weights = nullptr;

    return *this;
}


HotspotNucleus::~HotspotNucleus()
{
    safe_delete_hotspot_pos();
}


void HotspotNucleus::safe_delete_hotspot_pos()
{
    if (nullptr != m_hotspot_pos)
    {
        delete[] m_hotspot_pos;
        m_hotspot_pos = nullptr;
    }
}


void HotspotNucleus::prepare_hotspot_pos()
{
    safe_delete_hotspot_pos();

    m_hotspot_pos = new(std::nothrow) HotspotPos [m_num_hotspots_total];
    if (nullptr == m_hotspot_pos)
        exit(32);
}


void HotspotNucleus::safe_delete_hotspot_weights()
{
    if (nullptr != m_hotspot_weights)
    {
        delete[] m_hotspot_weights;
        m_hotspot_weights = nullptr;
    }
}


void HotspotNucleus::prepare_hotspot_weights()
{
    safe_delete_hotspot_weights();

    m_hotspot_weights = new(std::nothrow) double [m_num_hotspots_total];
    if (nullptr == m_hotspot_weights)
        exit(32);
}


void HotspotNucleus::sample_single_hotspot_pos (HotspotPos* hotspot_pos)
{
    hotspot_pos->x = m_rand_gaussian();
    hotspot_pos->y = m_rand_gaussian();
}


void HotspotNucleus::sample_hotspots_single_nucleon (uint nucleon_num)
{
    HotspotPos center_of_mass = {0.0, 0.0};

    for (uint i=0; i<m_num_hotspots_per_nucleon; i++)
    {
        HotspotPos* hotspot_pos = (HotspotPos*)get_hotspot_pos(nucleon_num, i);
        sample_single_hotspot_pos(hotspot_pos);

        center_of_mass.x += hotspot_pos->x;
        center_of_mass.y += hotspot_pos->y;
    }

    double inverse_divisor = 1.0/double(m_num_hotspots_per_nucleon);
    center_of_mass.x *= inverse_divisor;
    center_of_mass.y *= inverse_divisor;

    for (uint i=0; i<m_num_hotspots_per_nucleon; ++i)
    {
        HotspotPos* hotspot_pos = (HotspotPos*)get_hotspot_pos(nucleon_num, i);

        hotspot_pos->x -= center_of_mass.x;
        hotspot_pos->y -= center_of_mass.y;

        hotspot_pos->x += m_nucleon_pos[nucleon_num].x;
        hotspot_pos->y += m_nucleon_pos[nucleon_num].y;
    }
}


bool HotspotNucleus::fits_hotspot_distribution (double r_sqr)
{
    return m_rand() < NormedSamplingDistributions::gaussian(r_sqr, m_nucleon_size*m_nucleon_size);
}