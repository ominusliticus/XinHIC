#!/bin/python3

from numpy import pi, inf
from numpy import exp, sinh, cosh, log
from numpy import sqrt, fabs
from scipy.integrate import quad

# Angle averages of non-relativistic Maxwell-Boltzmann distirbution
def AngleAverageF(p, k, params):
    """ Eq. (A5) from paper """
    beta = 1 / params.T
    m = params.pion_mass
    r = m / params.D_mass
    
    x = beta * r * p * k / m
    return exp(-beta * (k**2 + (r * p)**2) / (2 * m)) * (sinh(x) / x)

def AngleAverageG(p, k, params):
    pass

# Real and Imaginary parts of D-mesons self energies
def ReDmesonSelfEnergy(p, params, pion_density, flag):
    """ 
    Real part of D-meson self energy 
    """
    m = params.pion_mass
    M = params.D_mass
    r = m / M
    Mt = M + m
    
    if flag == '0':
        d_1 = params.delta_zeroplus
        d_2 = params.delta_zerozero
    else: # flag must then be '+'
        d_1 = params.delta_plusminus
        d_2 = params.delta_pluszero

    def integrand(k, delta, params):
        """ 
        Integrand for principal value integral: Eq. (A6) 
        The small difference limit was obtained by leading 
        order expansion of Eq. (A6). I think this should be good
        enough
        """
        m = params.pion_mass
        M = params.D_mass
        r = m / M
        Mt = M + m
        
        # Bug for when delta is negative
        if delta < 0:
            delta = -delta
        K = sqrt(2 * r * Mt * delta)
        norm = 0.5 * K * ((2 * pi) / (m * params.T)) ** (1.5) / pi**2
        
        if fabs(k - K) < 5e-2:
            beta = 1 / params.T
            m = params.pion_mass
            r = m / params.D_mass
            
            x = beta * r * p * K / m
            val = exp(-beta * (K**2 + (r * p)**2) / (2 * m))
            val = 0.5 * val * (m * x * cosh(x) + (m * K**2 * beta) * sinh(x)) / (m * x)
        else:
            val = (k**2 * AngleAverageF(p, k, params) - K**2 * AngleAverageF(p, K, params))
            val = val / (k**2 - K**2)
        
        return norm * val

    mu_pi = m * M / (M + m)
    g_fpi2 = params.g**2 / params.f_pi**2
    K_1 = 2 * r * Mt * d_1
    K_2 = 2 * r * Mt * d_2 
    return - (mu_pi / m) * g_fpi2 * pion_density * (
            1.5 + (quad(integrand, 0, K_1, args=(d_1,params))[0] + \
                    quad(integrand, K_1, inf, args=(d_1,params))[0]) + \
                        0.5 * (quad(integrand, 0, K_2, args=(d_2,params))[0] + \
                            quad(integrand, K_2, inf, args=(d_2,params))[0]))   

def ImDmesonSelfEnergy(p, params, pion_density, flag):
    """
    Imaginary part of the D-meson self enery
    """
    m = params.pion_mass
    M = params.D_mass
    r = m / M
    Mt = M + m
    
    if flag == '0':
        d_1 = params.delta_zeroplus
        d_2 = params.delta_zerozero
    else: # flag must be then be '+'
        d_1 = -params.delta_plusminus
        d_2 = params.delta_pluszero
    
    mu_pi = m * M / (M + m)
    g_fpi2 = params.g**2 / params.f_pi**2 
    K = 2 * r * Mt 
    val = (mu_pi / (4 * pi * m)) * g_fpi2 * pion_density * (2 * K * pi / (params.T * m)) ** (1.5)
    val = val * (d_1**(1.5) * AngleAverageF(p, sqrt(K * d_1), params) + \
            0.5 * d_2**(1.5) * AngleAverageF(p, sqrt(K * d_2), params))
    return val


# Real and Imaginary part D*-mesons self energies
def ReDtmesonSelfEnergy(p, params):
    pass

def ImDtmesonSelfEnergy(p, parms):
    pass
