//
// Author: Kevin Ingles
// File: main.cpp

#include <fmt/core.h>
#include <cmath>
#include <fstream>
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
        size_t steps = static_cast<size_t>(std::floor((pmax - pmin) / 0.1));
        for (size_t i = 1; i < steps; i++)
        {
            double p = pmin + static_cast<double>(i) * deltap;
            fout << p << " ";
            for (auto temperature : temperatures)
                fout << ReDSelfEnergy(p, temperature, params.d0p, params.d00) << " " << ImDSelfEnergy(p, temperature, params.d0p, params.d00) << " ";
            fout << '\n';
        }
        fout.close();
    }
    return 0;
}
