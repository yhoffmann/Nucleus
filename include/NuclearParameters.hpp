#pragma once


#include <stdlib.h>


namespace NuclearParameters
{
    double mean_bulk_radius (uint atomic_num) // TODO find sources for these values
    {
        switch (atomic_num)
        {
            case 1:
                return 1.0; // fm

            case 16: // O
                return 2.608; // fm
            
            case 40: // Ar
                return 3.53;

            case 63: // Cu
                return 4.2;

            case 129: // Xe
                return 5.42;

            case 197: // Au
                return 6.38;

            case 208: // Pb
                return 6.62;

            default:
                exit(30);
        }
    }


    double mean_surface_diffusiveness (uint atomic_num)
    {
        switch (atomic_num)
        {
            case 1:
                return 0.01; // fm

            case 16: // O
                return 0.513; // fm
            
            case 40: // Ar
                return 0.542;

            case 63: // Cu
                return 0.513;

            case 129: // Xe
                return 0.57;

            case 197: // Au
                return 0.535;

            case 208: // Pb
                return 0.546;

            default:
                exit(30);
        }
    }
}