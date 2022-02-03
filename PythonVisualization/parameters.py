#!/bin/pyhton3

# class for containing all the necessary parameters in 
# our calculations

class Parameters:
    def __init__(self):
        # coupling and decay constant
        self.g = (0.329) ** 0.5
        self.f_pi = 130.5
        
        # PDG masses
        self.m_pion_plus = 139.57039
        self.m_pion_zero = 134.9758
        self.m_D_plus = 1869.58
        self.m_D_zero = 1864.83
        self.m_Dt_plus = 2010.26
        self.m_Dt_zero = 2006.85
        self.m_X = 3871.82

        # XEFT masses
        self.pion_mass = self.m_pion_zero
        self.D_mass = self.m_D_zero

        # Rest energies
        self.ep_pion_plus = self.m_pion_plus - self.pion_mass
        self.ep_pion_zero = self.m_pion_zero - self.pion_mass
        self.ep_D_plus = self.m_D_plus - self.D_mass
        self.ep_D_zero = self.m_D_zero - self.D_mass
        self.ep_Dt_plus = self.m_Dt_plus - (self.D_mass + self.pion_mass)
        self.ep_Dt_zero = self.m_Dt_zero - (self.D_mass + self.pion_mass)

        # difference in rest energies
        self.delta_pluszero = self.ep_Dt_plus - self.ep_D_plus - self.ep_pion_zero
        self.delta_zeroplus = self.ep_Dt_plus - self.ep_D_zero - self.ep_pion_plus
        self.delta_plusminus = self.ep_Dt_zero - self.ep_D_plus - self.ep_pion_plus
        self.delta_zerozero = self.ep_Dt_zero - self.ep_D_zero - self.ep_pion_zero

        # X binding energy
        self.ep_X = self.m_X - self.m_Dt_zero - self.m_D_zero
        self.gamma_X = (2 * self.D_mass * \
                (self.D_mass + self.pion_mass) / (2 * self.D_mass + self.pion_mass)\
                * self.ep_X) ** 0.5

        self.T = 0
