#include "../include/HotspotNucleus.hpp"


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

    m_hotspot_pos = new double [2*m_atomic_num*m_num_hotspots_per_nucleon];
    if (m_hotspot_pos == nullptr)
        exit(32);
}


void HotspotNucleus::sample_single_hotspot_pos (double* hotspot_pos)
{
    while (true)
    {
        hotspot_pos[0] = m_rand_gaussian();
        hotspot_pos[1] = m_rand_gaussian();

        double r_sqr = hotspot_pos[0]*hotspot_pos[0]+hotspot_pos[1]*hotspot_pos[1];

        if (m_rand() < NormedSamplingDistributions::gaussian(r_sqr, m_nucleon_size))
            return;
    }
}


void HotspotNucleus::sample_hotspots_single_nucleon (uint nucleon_num)
{
    double center_of_mass[2] = {0.0, 0.0};

    for (uint i=0; i<m_num_hotspots_per_nucleon; i++)
    {
        double* hotspot_pos = m_hotspot_pos+2*nucleon_num*m_num_hotspots_per_nucleon+2*i;
        sample_single_hotspot_pos(hotspot_pos);

        center_of_mass[0] += hotspot_pos[0];
        center_of_mass[1] += hotspot_pos[1];
    }

    center_of_mass[0] /= double(m_num_hotspots_per_nucleon);
    center_of_mass[1] /= double(m_num_hotspots_per_nucleon);

    for (uint i=0; i<m_num_hotspots_per_nucleon; i++)
    {
        double* hotspot_pos = m_hotspot_pos+2*nucleon_num*m_num_hotspots_per_nucleon+2*i;

        hotspot_pos[0] -= center_of_mass[0];
        hotspot_pos[1] -= center_of_mass[1];

        hotspot_pos[0] += m_nucleon_pos[3*nucleon_num];
        hotspot_pos[1] += m_nucleon_pos[3*nucleon_num+1];
    }
}


bool HotspotNucleus::fits_distribution (double r_sqr)
{
    return m_rand() < NormedSamplingDistributions::gaussian(r_sqr, m_nucleon_size);
}


void HotspotNucleus::set_hotspot_size (double hotspot_size)
{
    m_hotspot_size = hotspot_size;
}


void HotspotNucleus::sample_hotspots()
{
    for (uint n=0; n<m_atomic_num; n++)
    {
        sample_hotspots_single_nucleon(n);
    }
}


uint HotspotNucleus::get_num_hotspots_per_nucleon() const
{
    return m_num_hotspots_per_nucleon;
}


const double* HotspotNucleus::get_hotspot_pos (uint nucleon_num, uint hotspot_num) const
{
    return m_hotspot_pos+2*nucleon_num*m_num_hotspots_per_nucleon+2*hotspot_num;
}


HotspotNucleus::HotspotNucleus (uint atomic_num, uint num_hotspots_per_nucleon, std::mt19937& rng, SamplingDistribution sampling_distribution)
    : Nucleus(atomic_num, rng, sampling_distribution), m_num_hotspots_per_nucleon(num_hotspots_per_nucleon)
{
    prepare_hotspot_pos();
    sample_hotspots();

    m_dist_gaussian = std::normal_distribution<double>(0.0, m_nucleon_size);
}


HotspotNucleus::~HotspotNucleus()
{
    safe_delete_hotspot_pos();
}