//
// Author: Kevin Ingles
// File: MeatAndBones.hpp
// Description: This code contains all the functions that are needed 
//              to calculate evolution of X in warm hadron gas

#ifndef MEAT_AND_BONES_X_IN_HIC_HPP
#define MEAT_AND_BONES_X_IN_HIC_HPP   

namespace xinhic
{
    struct Parameters
    {
        // Note: All unitful quantities should be in terms MeV
        double g2   = 0.329;
        double fpi  = 130.5;

        // PDG mass
        double mplus    = 139.57039;
        double mzero    = 134.9768;
        double Mplus    = 1869.58;
        double Mzero    = 1864.83;
        double Mtplus   = 2010.26;
        double Mtzero   = 2006.85;
        double MT_true  = 3874.68;
        double MX_true  = 3871.82;

        // Galilean-invariant kinetic masses
        double m    = mzero;
        double M    = Mzero;
        double Mt   = M + m;
        double MX   = 2.0 * M + m;

        // Rest energies
        double epi0 = mzero - m;
        double epip = mplus - m;
        double e0   = Mzero - M;
        double ep   = Mplus - M;
        double et0  = Mtzero - Mt;
        double etp  = Mtplus - Mt;
        double eX   = MX_true - MX;
        double eT   = MT_true - 2.0 * Mzero - mplus;

        // Rest differences
        double d00 = et0 - e0 - epi0;
        double d0p = etp - e0 - epip;
        double dp0 = etp - ep - epi0;
        double dpm = et0 - ep - epip;
    };

    double ReDSelfEnergy(double D_energy, double temperature, double delta1, double delta2, double pion_density);
    double ImDSelfEnergy(double D_energy, double temperature, double delta1, double delta2, double pion_density);
    double ReDtSelfEnergyParallel(double Dt_energy, double temperature, double delta1, double delta2, double pion_density);  
    double ReDtSelfEnergyTransverse(double Dt_energy, double temperature, double delta1, double delta2, double pion_density);
    double ImDtSelfEnergyParallel(double Dt_energy, double temperature, double delta1, double delta2, double pion_density);
    double ImDtSelfEnergyTransverse(double Dt_energy, double temperature, double delta1, double delta2, double pion_density);
}

#endif
