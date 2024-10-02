/*
sources:

[1] 
    Improved version of the PHOBOS Glauber Monte Carlo
    Loizides, C. et al.
    SoftwareX, Volume 1, 13 - 18, 2015

[2]
    Improved Monte Carlo Glauber predictions at present and future nuclear colliders
    Constantin Loizides, Jason Kamin, and David d'Enterria
    Phys. Rev. C 97, 054910, 2018
*/

#pragma once


#include <stdlib.h>


const char* ERROR_NO_DATA = "No data for the nucleon number. Exiting!";


struct NucleusParameters
{
    double mean_bulk_radius; // R
    double mean_surface_diffusiveness; // a
    double deviation_from_spherical; // w
};

namespace NuclearParameters
{
    // from source [1] and [2]
    const NucleusParameters params[] =
    { // R,     a,      w               Symbol  A           index
        {1.0,   0.0001, 0.0},       //  p       1           0       // surface diffusiveness small to speed up the "sampling", the value itself does not mean anything
        {2.608, 0.513,  -0.51},     //  O       16          1
        {3.34,  0.58,   -0.233},    //  Si      28          2
        {2.54,  2.191,  0.16},      //  S       32          3
        {3.766, 0.586,  -0.161},    //  Ca      40          4
        {4.309, 0.517,  -0.1308},   //  Ni      58          5
        {4.2,   0.596,  0.0},       //  Cu      62,63       6
        {5.36,  0.59,   0.0},       //  Xe      129         7
        {6.58,  0.48,   0.0},       //  W       186         8
        {6.38,  0.535,  0.0},       //  Au      197         9
        {6.62,  0.546,  0.0},       //  Pb      207,208     10
    };

    inline const NucleusParameters& get (uint atomic_num)
    {
        switch (atomic_num)
        {
            case 1:
                return params[0];

            case 16:
                return params[1];

            case 28:
                return params[2];

            case 32:
                return params[3];

            case 40:
                return params[4];

            case 58:
                return params[5];

            case 62:
            case 63:
                return params[6];

            case 129:
                return params[7];

            case 186:
                return params[8];

            case 197:
                return params[9];

            case 207:
            case 208:
                return params[10];

            default:
                std::cerr << ERROR_NO_DATA << std::endl;
                exit(30);
        }
    }
}