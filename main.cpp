//
// Author: Kevin Ingles
// File: main.cpp

#include <fmt/core.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "MeatAndBones.hpp"

using namespace xinhic;

static Parameters params{};

int main()
{
    double pmin = 0.0;
    double pmax = 120.0;
    double deltap = 0.1;

    {
        double temperatures[5] = {150.0, 125.0, 100.0, 75.0, 50.0}; // MeV
        auto fout = std::fstream("trial_file.dat", std::fstream::out);
        size_t steps = static_cast<size_t>(std::floor((pmax - pmin) / deltap));
        for (size_t i = 1; i < steps; i++)
        {
            double p = pmin + static_cast<double>(i) * deltap;
            fout << p << " ";
            for (auto temperature : temperatures)
                fout << ReDSelfEnergy(p, temperature, params.d0p, params.d00, 2.3e5) << " " << ImDSelfEnergy(p, temperature, params.d0p, params.d00, 2.3e5) << " ";
            fout << '\n';
        }
        fout.close();
    }
    const double PI = 4.0 * std::atan(1.0);
    double G_pi2 = params.g2 / (4.0 * PI * params.m * params.fpi * params.fpi);
    double mu_piT = params.m * params.MX / (2.0 * params.Mt);
    double mu_pi = params.m * params.M / params.Mt;
    std::cout << " G_pi^2 mu_piT mu_pi m = " << G_pi2 * mu_piT * mu_pi * params.m << '\n';
    return 0;
}
