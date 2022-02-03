//
// Author: Kevin Ingles
// File: main.cpp

#include <cmath>
#include <iostream>
#include <fmt/core.h>
#include <fstream>
#include "MeatAndBones.hpp"
#include <utility>

template<class OStream, typename...Args>
void Print(OStream& os, Args...args)
{
    ((os << std::forward<Args>(args) << " "),...);
}

using namespace xinhic;

static Parameters params{};
static constexpr char newl = '\n';

int main()
{
    double Emin = 0.0;
    double Emax = 115.0;
    double deltaE = 0.1;
    
    {
        double pion_density = 0.016 * std::pow(200.0, 3.0);
        double temperatures[] = { 115 }; // MeV
        auto fout_D = std::fstream("D-self-energies.dat", std::fstream::out);
        auto fout_Dt = std::fstream("Dt-self-energies.dat", std::fstream::out);
        size_t steps = static_cast<size_t>(std::floor((Emax - Emin) / deltaE));
        for (size_t i = 1; i < steps; i++)
        {
            double E = Emin + static_cast<double>(i) * deltaE;
            Print(fout_D, E);
            Print(fout_Dt, E);
            for (auto temperature : temperatures)
            {
                Print(fout_D, ReDSelfEnergy(E, temperature, params.dp0, params.d00, pion_density), ImDSelfEnergy(E, temperature, params.dp0, params.d00, pion_density));
                Print(fout_D, ReDSelfEnergy(E, temperature, params.dpm, params.d0p, pion_density), ImDSelfEnergy(E, temperature, params.dpm, params.d0p, pion_density));

                Print(fout_Dt, ReDtSelfEnergyParallel(E, temperature, params.dpm, params.d00, pion_density), ReDtSelfEnergyTransverse(E, temperature, params.dpm, params.d00, pion_density));
                Print(fout_Dt, ImDtSelfEnergyParallel(E, temperature, params.dpm, params.d00, pion_density), ImDtSelfEnergyTransverse(E, temperature, params.dpm, params.d00, pion_density));
                Print(fout_Dt, ReDtSelfEnergyParallel(E, temperature, params.dp0, params.d0p, pion_density), ReDtSelfEnergyTransverse(E, temperature, params.dp0, params.d0p, pion_density));
                Print(fout_Dt, ReDtSelfEnergyParallel(E, temperature, params.dp0, params.d0p, pion_density), ReDtSelfEnergyTransverse(E, temperature, params.dp0, params.d0p, pion_density));
            }
            if (i < steps - 1)
            {
                Print(fout_D, newl);
                Print(fout_Dt, newl);
            }
        }
        fout_D.close();
        fout_Dt.close();
    }
    const double PI = 4.0 * std::atan(1.0);
    double G_pi2 = params.g2 / (4.0 * PI * params.m * params.fpi * params.fpi);
    double mu_piT = params.m * params.MX / (2.0 * params.Mt);
    double mu_pi = params.m * params.M / params.Mt;
    std::cout << " G_pi^2 mu_piT mu_pi m = " << G_pi2 * mu_piT * mu_pi * params.m << '\n';
    return 0;
}
