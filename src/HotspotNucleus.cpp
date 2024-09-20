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


void HotspotNucleus::sample_hotspot_pos()
{
    for (uint n=0; n<m_atomic_num; n++)
    {
        sample_hotspots_single_nucleon(n);
    }
}


void HotspotNucleus::sample_nucleon_pos()
{
    Nucleus::sample_nucleon_pos();

    sample_hotspot_pos();
}


double HotspotNucleus::get_hotspot_thickness (double x, double y) const
{
    double thickness = 0.0;

    double inverse_r_sqr_divisor = 1.0/(2.0*m_hotspot_size*m_hotspot_size);
    for (uint i=0, i_max=m_atomic_num*m_num_hotspots_per_nucleon; i<i_max; ++i)
    {
        double delta_x = x-m_hotspot_pos[i].x;
        double delta_y = y-m_hotspot_pos[i].y;

        double r_sqr = delta_x*delta_x + delta_y*delta_y;

        thickness += exp( -r_sqr*inverse_r_sqr_divisor );
    }

    return thickness/double(m_num_hotspots_per_nucleon)*inverse_r_sqr_divisor/M_PI;
}


uint HotspotNucleus::get_num_hotspots_per_nucleon() const
{
    return m_num_hotspots_per_nucleon;
}


const HotspotPos* HotspotNucleus::get_hotspot_pos (uint nucleon_num, uint hotspot_num) const
{
    return m_hotspot_pos+nucleon_num*m_num_hotspots_per_nucleon+hotspot_num;
}


HotspotNucleus::HotspotNucleus (uint atomic_num, uint num_hotspots_per_nucleon, uint seed, SamplingDistribution sampling_distribution)
    : Nucleus(atomic_num, seed, sampling_distribution)
    , m_num_hotspots_per_nucleon(num_hotspots_per_nucleon)
{   
    prepare_hotspot_pos();
    sample_hotspot_pos();
}


HotspotNucleus::HotspotNucleus (const HotspotNucleus& other)
    : Nucleus(other)
    , m_num_hotspots_per_nucleon(other.m_num_hotspots_per_nucleon)
    , m_hotspot_size(other.m_hotspot_size)
{
    prepare_hotspot_pos();
    std::copy(other.m_hotspot_pos, other.m_hotspot_pos + m_atomic_num*m_num_hotspots_per_nucleon, m_hotspot_pos);
}


HotspotNucleus::HotspotNucleus (HotspotNucleus&& other)
    : Nucleus(std::move(other))
    , m_num_hotspots_per_nucleon(other.m_num_hotspots_per_nucleon)
    , m_hotspot_pos(other.m_hotspot_pos)
    , m_hotspot_size(other.m_hotspot_size)
{
    other.m_hotspot_pos = nullptr;
}


HotspotNucleus& HotspotNucleus::operator= (const HotspotNucleus& other)
{
    if (this == &other)
        return *this;

    Nucleus::operator=(other);
    m_num_hotspots_per_nucleon = other.m_num_hotspots_per_nucleon;
    m_hotspot_pos = other.m_hotspot_pos;
    m_hotspot_size = other.m_hotspot_size;

    prepare_hotspot_pos();
    std::copy(other.m_hotspot_pos, other.m_hotspot_pos + m_atomic_num*m_num_hotspots_per_nucleon, m_hotspot_pos);

    return *this;
}


HotspotNucleus& HotspotNucleus::operator= (HotspotNucleus&& other)
{
    if (this == &other)
        return *this;

    Nucleus::operator=(std::move(other));
    m_num_hotspots_per_nucleon = other.m_num_hotspots_per_nucleon;
    m_hotspot_pos = other.m_hotspot_pos;
    m_hotspot_size = other.m_hotspot_size;

    other.m_hotspot_pos = nullptr;

    return *this;
}


HotspotNucleus::~HotspotNucleus()
{
    safe_delete_hotspot_pos();
}


void HotspotNucleus::safe_delete_hotspot_pos()
{
    if (m_hotspot_pos)
    {
        delete[] m_hotspot_pos;
        m_hotspot_pos = nullptr;
    }
}


void HotspotNucleus::prepare_hotspot_pos()
{
    safe_delete_hotspot_pos();

    m_hotspot_pos = new(std::nothrow) HotspotPos [m_atomic_num*m_num_hotspots_per_nucleon];
    if (m_hotspot_pos == nullptr)
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