
// Author: Kevin Ingles
// File: MeatAndBones.cpp
// Description: Translation unit for MeatAndBones.hpp

#include "MeatAndBones.hpp"

#include <limits>
#include <cmath>
#include <iostream>

#include <fmt/core.h>

static constexpr bool DEBUG = false;

static constexpr double inf = std::numeric_limits<double>::infinity();

// Choice to generalize code to calculate number of points to use
static constexpr int NSUM48 = 24;

// constexpr keyword is to ensure there are not defined in multiple translation units
static constexpr double x48[NSUM48] = { 0.0323801709628694,
                    0.0970046992094627,
                    0.1612223560688917,
                    0.2247637903946891,
                    0.2873624873554556,
                    0.3487558862921608,
                    0.4086864819907167,
                    0.4669029047509584,
                    0.5231609747222330,
                    0.5772247260839727,
                    0.6288673967765136,
                    0.6778723796326639,
                    0.7240341309238146,
                    0.7671590325157404,
                    0.8070662040294426,
                    0.8435882616243935,
                    0.8765720202742479,
                    0.9058791367155696,
                    0.9313866907065543,
                    0.9529877031604309,
                    0.9705915925462473,
                    0.9841245837228269,
                    0.9935301722663508,
                    0.9987710072524261 };

static constexpr double w48[NSUM48] = { 0.0647376968126839,
                    0.0644661644359501,
                    0.0639242385846482,
                    0.0631141922862540,
                    0.0620394231598927,
                    0.0607044391658939,
                    0.0591148396983956,
                    0.0572772921004032,
                    0.0551995036999842,
                    0.0528901894851937,
                    0.0503590355538545,
                    0.0476166584924905,
                    0.0446745608566943,
                    0.0415450829434647,
                    0.0382413510658307,
                    0.0347772225647704,
                    0.0311672278327981,
                    0.0274265097083569,
                    0.0235707608393244,
                    0.0196161604573555,
                    0.0155793157229438,
                    0.0114772345792345,
                    0.0073275539012763,
                    0.0031533460523058 };


///////////////////////////////////////////////////////
//     Prototyping so the functions see each other   //
///////////////////////////////////////////////////////
template<typename Functor, typename...Args>
static double GausQuadAux(Functor&& func, double _low, double _high, double result, double tol, int depth, bool improper_top, Args&&... args);

template<typename Functor, typename...Args>
static double GausQuad(Functor&& func, double _low, double _high, double tol, int maxDepth, Args&&... args);

///////////////////////////////////////////////////////
//              Defining implementation              //
///////////////////////////////////////////////////////
template<typename Functor, typename...Args>
static double GausQuadAux(Functor&& func, double _low, double _high, double result, double tol, int depth, bool improper_top, Args&&... args)
{
    // Quick check to ensure that we should do calculation
    if (depth <= 0)
    {
        if (DEBUG)
        {
            std::cerr << "Integrate failed to converge.\n";
            std::cerr << "Either increase depth or analyze integrand.\n";
        }
        return result;
    }
    else if (std::isnan(result))
    {
        std::cerr << "NaN encountered in integration, please check integrand.\n";
        // exit(-1);
    }

    double high = _high;
    double low = _low;

    double middle = (high + low) / 2.0;
    double interval1_result = 0;
    double interval2_result = 0;

    double yneg, ypos;
    for (int i = 0; i < NSUM48; i++)
    {
        // Sum up areas using above weights and points
        yneg = ((middle - low) * (-x48[i]) + (middle + low)) / 2.0;
        ypos = ((middle - low) * (x48[i]) + (middle + low)) / 2.0;
        if (!improper_top)
            interval1_result += w48[i] * func(yneg, std::forward<Args>(args)...) + w48[i] * func(ypos, std::forward<Args>(args)...);
        else
            interval1_result += w48[i] * func(1 / yneg, std::forward<Args>(args)...) / pow(yneg, 2.0) + w48[i] * func(1 / ypos, std::forward<Args>(args)...) / pow(ypos, 2.0);

        // Sum up areas using above weights and points
        yneg = ((high - middle) * (-x48[i]) + (high + middle)) / 2.0;
        ypos = ((high - middle) * (x48[i]) + (high + middle)) / 2.0;
        if (!improper_top)
            interval2_result += w48[i] * func(yneg, std::forward<Args>(args)...) + w48[i] * func(ypos, std::forward<Args>(args)...);
        else
            interval2_result += w48[i] * func(1 / yneg, std::forward<Args>(args)...) / pow(yneg, 2.0) + w48[i] * func(1 / ypos, std::forward<Args>(args)...) / pow(ypos, 2.0);
    }
    interval1_result *= (middle - low) / 2.0;
    interval2_result *= (high - middle) / 2.0;

    double result2 = interval1_result + interval2_result;

    if (fabs(result - result2) / result < tol)
        return result;
    else
        return GausQuadAux(func, low, middle, interval1_result, tol, depth - 1, improper_top, std::forward<Args>(args)...) + GausQuadAux(func, middle, high, interval2_result, tol, depth - 1, improper_top, std::forward<Args>(args)...);
}

template<typename Functor, typename...Args>
static double GausQuad(Functor&& func, double _low, double _high, double tol, int maxDepth, Args&&... args)
{
    double result = 0;
    double yneg, ypos;
    double high = _high;
    double low = _low;

    if (high == low) return 0.0;

    bool improper_top = false;
    if (high == inf || low == -inf)
    {
        if (high == inf && low != -inf)
        {
            if (low == 0)
            {
                result += GausQuad(func, 0, 1, tol, maxDepth, std::forward<Args>(args)...);
                result += GausQuad(func, 1, high, tol, maxDepth, std::forward<Args>(args)...);
                return result;
            }
            else high = 1 / low;

            low = 0;
            improper_top = true;
        }
        else if (high != inf && low == -inf)
        {
            if (high == 0)
            {
                result += GausQuad(func, -1, 0, tol, maxDepth, std::forward<Args>(args)...);
                result += GausQuad(func, low, -1, tol, maxDepth, std::forward<Args>(args)...);
                return result;
            }
            else low = 1 / high;

            high = 0;
            improper_top = true;
        }
        else
        {
            result += GausQuad(func, low, -1, tol, maxDepth, std::forward<Args>(args)...);
            result += GausQuad(func, -1, 1, tol, maxDepth, std::forward<Args>(args)...);
            result += GausQuad(func, 1, high, tol, maxDepth, std::forward<Args>(args)...);
            return result;
        }
    }

    for (int i = 0; i < NSUM48; i++)
    {
        // Sum up areas using above weights and points
        yneg = ((high - low) * (-x48[i]) + (high + low)) / 2.0;
        ypos = ((high - low) * (x48[i]) + (high + low)) / 2.0;
        if (!improper_top)
        {
            result += w48[i] * func(yneg, std::forward<Args>(args)...) + w48[i] * func(ypos, std::forward<Args>(args)...);
        }
        else
        {
            result += w48[i] * func(1 / yneg, std::forward<Args>(args)...) / pow(yneg, 2.0) + w48[i] * func(1 / ypos, std::forward<Args>(args)...) / pow(ypos, 2.0);
        }
    }
    result *= (high - low) / 2.0;

    return GausQuadAux(func, low, high, result, tol, maxDepth, improper_top, std::forward<Args>(args)...);
}

static double PI = 4.0 * std::atan(1.0);

namespace xinhic
{    
    static Parameters params{};

    static double F(double p, double k, double temperature, double r)
    {
        double T = temperature;
        double exp_part = std::exp(-(k * k + r * r * p * p) / (2.0 * params.m * T));
        double sinhc_part = (k * p == 0) ? 1 : std::sinh(r * k * p / (params.m * T)) / (r * k * p / (params.m * T));
        double value = exp_part * sinhc_part;
        return exp_part > 0 ? value : 0;
    }

    static double G(double p, double k, double temperature)
    {
        double T = temperature;
        double T2 = T * T;
        double T3 = T2 * T;
        double r = params.m / params.Mt;
        double exp_part = std::exp(-(k * k + r * r * p * p) / (2.0 * params.m * T));
        double sinhc_part = std::sinh(r * k * p / (params.m * T)) / r * k * p / (params.m * T);
        double coshc_part = std::cosh(r * k * p / (params.m * T)) / r * k * p / (params.m * T);
        double derivative = -2.0 * T2  * coshc_part;
        derivative += 2.0 * params.m * T3 * sinhc_part / (r * k * p);
        derivative += r * k * p * T * sinhc_part / params.m;
        return exp_part * std::pow(params.m / (r * p * k), 2.0) * derivative;
    }

    double ReDSelfEnergy(double D_energy, double temperature, double delta1, double delta2, double pion_density)
    {
        // Interand for self energy
        double p = std::sqrt(2.0 * params.M * D_energy);
        double T = temperature;
        double r = params.m / params.M;
        auto integrand = [=](double k, double delta) -> double
        {
            double k2 = k * k;
            double K2 = 2 * r * params.Mt * delta;
            if (delta < 0)
                return k2 * F(p, k, T, r) / (k2 - K2);
            else
            {
                double K = std::sqrt(K2);
                if (k2 - K2 == 0)
                    return 0.5 * F(K, K, T, r);
                else 
                    return (k2 * F(p, k, T, r) - K2 * F(p, K, T, r)) / (k2 - K2);
            }
        };
        double K1 = (delta1 > 0) ? std::sqrt(2.0 * r * params.Mt * delta1) : inf;
        double K2 = (delta2 > 0) ? std::sqrt(2.0 * r * params.Mt * delta2) : inf;
        double I1 = GausQuad(integrand, 0, K1, 1e-6, 10, delta1) + GausQuad(integrand, K1, inf, 1e-6, 10, delta1);
        double I2 = GausQuad(integrand, 0, K2, 1e-6, 10, delta2) + GausQuad(integrand, K2, inf, 1e-6, 10, delta2);
        
        double C1 = r * params.Mt * delta1;
        double C2 = r * params.Mt * delta2;
        double dist_norm = std::pow(2.0 * PI / (params.m * T), 1.5) / (PI * PI);
        double prefactor = - pion_density * params.g2 * params.M / (params.Mt * params.fpi * params.fpi);
        double return_val = 2.0 * pion_density / (params.fpi * params.fpi);
        return_val += prefactor * (1.5 + dist_norm * (C1 * I1 + 0.5 * C2 * I2));
        return return_val;
    }

    double ImDSelfEnergy(double D_energy, double temperature, double delta1, double delta2, double pion_density)
    {
        double p = std::sqrt(2.0 * params.M * D_energy);
        double T = temperature;
        double r = params.m / params.M;
        double x2 = 2.0 * r * params.Mt; // really can't think of a good descriptionve name
        double x = std::sqrt(x2);
        
        double dist_norm = std::pow(2.0 * PI / (params.m * T), 1.5) / (PI * PI);
        double prefactor = pion_density * PI * params.g2 * params.M / (4.0 * params.Mt * params.fpi * params.fpi);
        double return_value = 0.0;
        return_value += (delta1 > 0) ? std::pow(x2 * delta1, 1.5) * F(p, x * std::sqrt(delta1), T, r) : 0;
        return_value += (delta2 > 0) ? std::pow(x2 * delta2, 1.5) * F(p, x * std::sqrt(delta2), T, r) : 0;
        
        return prefactor * dist_norm * return_value;
    }

    double ReDtSelfEnergyParallel(double Dt_energy, double temperature, double delta1, double delta2, double pion_density)
    {
        double q = std::sqrt(2.0 * (params.M + params.m) * Dt_energy);
        double T = temperature;
        double m = params.m;
        double M = params.M;
        double Mt = M + m;
        double r = params.m / Mt;

        auto integrandF = [=](double k, double delta) -> double
        {
            double k2 = k * k;
            double K2 = 2 * (params.m * params.M) / Mt * delta;
            if (delta < 0)
                return k2 * F(q, k, T, r) / (k2 - K2);
            else
            {
                double K = std::sqrt(K2);
                if (k2 - K2 == 0)
                    return 0.5 * F(K, K, T, r);
                else 
                    return (k2 * F(q, k, T, r) - K2 * F(q, K, T, r)) / (k2 - K2);
            }
        };

        auto integrandG = [=](double k, double delta) -> double
        {
            double k2 = k * k;
            double K2 = 2 * (params.m * params.M) / Mt * delta;
            if (delta < 0)
                return k2 * G(q, k, T) / (k2 - K2);
            else
            {
                double K = std::sqrt(K2);
                if (k2 - K2 == 0)
                    return 0.5 * G(K, K, T);
                else 
                    return (k2 * G(q, k, T) - K2 * G(q, K, T)) / (k2 - K2);
            }
        };

        auto integrand_log = [=](double k)
        {
            if (k == 0) return 0.0;
            if (k == r * q) return 0.0;
            double r2 = r * r;
            double r3 = r2 * r;
            double k2 = k * k;
            double q2 = q * q;
            double q3 = q2 * q;
            double particle_dist = std::exp(- k2 / (2.0 * m * T));
            double frac_1 = std::pow(k2 - r2 * q2, 2.0) / (16.0 * r3 * k * q3);
            double frac_2 = std::pow((k + r * q) / (k - r * q), 2.0);
            return k2 * particle_dist * frac_1 * std::log(frac_2);
        };

        double K1 = (delta1 > 0) ? std::sqrt(2.0 * r * params.Mt * delta1) : inf;
        double K2 = (delta2 > 0) ? std::sqrt(2.0 * r * params.Mt * delta2) : inf;
        double IF1 = GausQuad(integrandF, 0, K1, 1e-6, 10, delta1) + GausQuad(integrandF, K1, inf, 1e-6, 10, delta1);
        double IF2 = GausQuad(integrandF, 0, K2, 1e-6, 10, delta2) + GausQuad(integrandF, K2, inf, 1e-6, 10, delta2);
        double IG1 = GausQuad(integrandG, 0, K1, 1e-6, 10, delta1) + GausQuad(integrandG, K1, inf, 1e-6, 10, delta1);
        double IG2 = GausQuad(integrandG, 0, K2, 1e-6, 10, delta2) + GausQuad(integrandG, K2, inf, 1e-6, 10, delta2);

        double I_log = GausQuad(integrand_log, 0, r * q, 1e-6, 10) + GausQuad(integrand_log, r * q, inf, 1e-6, 10);

        double C1 = r * params.Mt * delta1;
        double C2 = r * params.Mt * delta2;
        double dist_norm = std::pow(2.0 * PI / (params.m * T), 1.5) / (PI * PI);
        double prefactor = params.g2 / (2.0 * params.fpi * params.fpi) * std::pow(M / Mt, 3.0);
        
        double A1 = 0.5 * C1 * dist_norm * (IF1 - IG1);
        double A2 = 0.5 * C2 * dist_norm * (IF2 - IG2);

        double return_val =  0.25 * (1.0 - 6.0 * params.m * T / (r * r * q * q)) * pion_density;
        return_val += 0.25 * dist_norm * I_log;
        return_val += 0.5 * (A1 + A2);
        return prefactor * return_val + 2.0 * pion_density / (params.fpi * params.fpi);                
    }

    double ReDtSelfEnergyTransverse(double Dt_energy, double temperature, double delta1, double delta2, double pion_density)
    {
        return -0.5 * (ReDtSelfEnergyParallel(Dt_energy, temperature, delta1, delta2, pion_density)
                - 2.0 * pion_density / (params.fpi * params.fpi));
    }
    
    double ImDtSelfEnergyParallel(double Dt_energy, double temperature, double delta1, double delta2, double pion_density)
    {
        double m = params.m;
        double M = params.M;
        double Mt = M + m;
        double r = m / Mt;
        double mu_pi = m * M / Mt;
        double T = temperature;

        double q = std::sqrt(2.0 * Mt * Dt_energy);

        double prefactor = params.g2 / (2.0 * params.fpi * params.fpi) * std::pow(M / Mt, 3.0);
        double dist_norm = std::pow(2.0 * PI / (m * T), 1.5) / PI;

        double C1 = delta1 > 0 ? std::sqrt(2.0 * mu_pi * delta1) : 0.0;
        double C2 = delta2 > 0 ? std::sqrt(2.0 * mu_pi * delta2) : 0.0;
        double return_val = std::pow(C1, 3.0) * (F(q, C1, T, r) - G(q, C1, T));
        return_val += std::pow(C2, 3.0) * (F(q, C2, T, r) - G(q, C2, T));        
        return prefactor * dist_norm * return_val;
    }
    
    double ImDtSelfEnergyTransverse(double Dt_energy, double temperature, double delta1, double delta2, double pion_density)
    {
        return -0.5 * ImDtSelfEnergyParallel(Dt_energy, temperature, delta1, delta2, pion_density);
    }
}
