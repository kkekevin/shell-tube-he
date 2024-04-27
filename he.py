import numpy as np

def tube_side_coefficient(d_i, m_dot, c_p, mu, k, l):
    # Calculate the Reynolds number (Re)
    Re = (4 * m_dot) / (np.pi * d_i * mu)
    
    # Calculate the Prandtl number (Pr) if not provided
    pr = c_p * mu / k

    # Calculate the Nusselt number using Dittus-Boelter equation (for turbulent flow)
    Nu = 0.023 * Re**0.8 * pr**0.4
    
    # Calculate the tube-side heat transfer coefficient (h_tube)
    h_tube = Nu * k / d_i
    
    return h_tube

def shell_side_coefficient(d_o, baffle_cut, pitch, m_dos, c_p, mu, k, l, D_shell):
    # Placeholder for actual shell-side calculation. Simplified version:
    # Calculate Reynolds number using an estimated cross-flow area
    # (assumes 25% of total shell area cut away by baffles)
    A_cross = np.pi * (D_shell**2 / 4) * (1 - baffle_cut)
    V_shell = m_dos / (A_cross * 780)  # assuming shell-side fluid density (water) is 780 kg/m3
    Re_shell = d_o * V_shell / mu
    pr = c_p * mu / k # not provided
    
    # Use Colburn and friction factor relations for idealized shell-side flow
    jH = 0.5 * Re_shell**(-0.2)  # Colburn analogy for moderate Reynolds numbers
    h_shell = jH * k / d_o * Re_shell * pr**(1/3)
    
    return h_shell

def overall_heat_transfer_coefficient(h_tube, h_shell, d_o, d_i, k_wall):
    # Wall resistance and fouling factors (simplified, no fouling)
    R_wall = np.log(d_o / d_i) / (2 * np.pi * k_wall)
    
    # Calculate overall heat transfer coefficient (U)
    U = 1 / ((1 / h_tube) + R_wall + (1 / h_shell) * (d_o / d_i))
    
    return U

def tlm (tempIT, tempOT, tempIS, tempOS) :
    return (tempIS - tempOT) - (tempOS - tempIT) / np.log((tempIS - tempOT) / (tempOS - tempIT))

def required_area(Q, U, deltaT_lm):
    # Calculate required heat transfer area
    A = Q / (U * deltaT_lm)
    return A

# tube side parameters
m_dot = 22.22  # mass flow rate in kg/s
c_p = 2200.8   # specific heat capacity in J/kgK (ketone)
mu = 0.4131  # dynamic viscosity in Pa.s
k = 0.6  # thermal conductivity in W/mK
tempInTube = 283.15
tempOutTube = 303.15
#pr = 6.0  # Prandtl number (for water)
l = 4    # tube length in meters
d_i = 0.015  # inner diameter of the tube in meters
d_o = 0.019  # outer diameter of the tube in meters
k_wall = 15  # thermal conductivity of tube wall in W/mK
# shell side parameters
m_dos = 27.77
c_pShell = 1900
muShell = 0.00152
tempInShell = 373.15
tempOutShell = 363.15
D_shell = 0.5  # shell diameter in meters
baffle_cut = 0.25  # baffle cut as a fraction of the shell diameter
Q = 50000  # total heat transfer in Watts
deltaT_lm = tlm(tempInTube, tempOutTube, tempInShell, tempOutShell)  # log mean temperature difference in K

# Calculate coefficients
h_tube = tube_side_coefficient(d_i, m_dot, c_p, mu, k, l)
h_shell = shell_side_coefficient(d_o, baffle_cut, d_o*1.25, m_dos, c_pShell, muShell, k, l, D_shell)
U = overall_heat_transfer_coefficient(h_tube, h_shell, d_o, d_i, k_wall)
A = required_area(Q, U, deltaT_lm)

print(f"Tube-side heat transfer coefficient: {h_tube:.2f} W/m²K")
print(f"Shell-side heat transfer coefficient: {h_shell:.2f} W/m²K")
print(f"Overall heat transfer coefficient: {U:.2f} W/m²K")
print(f"Required heat transfer area: {A:.2f} m²")