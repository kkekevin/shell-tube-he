import numpy as np
import math

def numOfTubes (shell_diameter, tube_outer_diameter, layout) :
    """
    Calculate the estimated number of tubes that can be fitted in a given shell diameter based on Kern's method.

    Parameters:
    - shell_diameter (float): Inner diameter of the shell (m).
    - tube_outer_diameter (float): Outer diameter of the tubes (m).
    - tube_pitch (float): Center-to-center distance between adjacent tubes (m).
    - layout (str): Tube layout pattern, either 'triangle' or 'square'.

    Returns:
    int: Estimated number of tubes that can fit within the shell.
    """
    if layout not in ['t', 's']:
        raise ValueError("Layout must be either 'triangle' or 'square'.")

    # Calculate the number of tube rows that can fit along the shell diameter
    # The effective shell diameter is reduced by one tube diameter to account for clearance
    effective_diameter = shell_diameter - tube_outer_diameter

    if layout == 't':
        tube_pitch = 1.25 * tube_outer_diameter
        # For triangular pitch, calculate side length of an equilateral triangle
        tube_rows = math.floor((2 / math.sqrt(3)) * (effective_diameter / tube_pitch))
        # The number of tubes in each row follows 1, 2, 3, ..., n pattern
        num_tubes = (tube_rows * (tube_rows + 1)) // 2
    elif layout == 's':
        tube_pitch = 1.5 * tube_outer_diameter
        # For square pitch, straightforward division
        tube_rows = math.floor(effective_diameter / tube_pitch)
        # Number of tubes per row is equal to the number of rows
        num_tubes = tube_rows ** 2

    return num_tubes

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

def shell_side_coefficient(d_o, baffle_cut, pitch, m_dos, c_p, mu, k, l, sg, D_shell):
    # Placeholder for actual shell-side calculation. Simplified version:
    # Calculate Reynolds number using an estimated cross-flow area
    # (assumes 25% of total shell area cut away by baffles)
    A_cross = np.pi * (D_shell**2 / 4) * (1 - baffle_cut)
    V_shell = m_dos / (A_cross * sg)  # assuming shell-side fluid density (ketone) is 780 kg/m3
    Re_shell = d_o * V_shell / mu
    pr = c_p * mu / k # not provided

    print(f'aaa: {V_shell: .2f}')
    print(Re_shell)
    
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

# log mean temperature for counter-flow
def tlm (tempInCold, tempOutCold, tempInHot, tempOutHot) :
    return ((tempInHot - tempOutCold) - (tempOutHot - tempInCold)) / np.log((tempInHot - tempOutCold) / (tempOutHot - tempInCold))

def required_area(Q, U, deltaT_lm):
    # Calculate required heat transfer area
    A = Q / (U * deltaT_lm)
    return A

def estimate_tube_number(heat_transfer_area, tube_outer_diameter, length, n):
    """
    Estimate the number of tubes required to achieve a specified heat transfer area without a predetermined shell diameter.
    
    Parameters:
    - heat_transfer_area (float): The total required heat transfer area (m^2).
    - tube_outer_diameter (float): The outer diameter of the tubes (m).
    - tube_pitch (float): The center-to-center distance between tubes (m).
    - layout (str): The layout of the tubes ('triangle' or 'square').
    
    Returns:
    int: The estimated total number of tubes.
    """
    # Calculate the surface area of one tube (assuming a unit length for simplicity)
    tube_surface_area_per_meter = math.pi * tube_outer_diameter  # per meter tube length

    # Test four length of tube
    for l in length :
        tube_surface_area = tube_surface_area_per_meter * l
        # Calculate the number of tubes needed to achieve the desired total heat transfer area
        total = math.ceil(heat_transfer_area / tube_surface_area)
        if (total <= n) :
            return [l, total]
    

# tube side parameters (therminol)
m_dot = 27.77  # mass flow rate in kg/s
c_p = 1900   # specific heat capacity in J/kgK (therminol)
mu = 0.00152  # dynamic viscosity in Pa.s
sg = 944 #density of the tube's fluid in kg/m³
k = 0.123  # thermal conductivity in W/mK (therminol60)
tempInTube = 373.15
tempOutTube = 370
#pr = 6.0  # Prandtl number (for water)
l = [1, 2, 3, 4]    # tube length in meters available
d_i = 0.015  # inner diameter of the tube in meters
d_o = 0.019  # outer diameter of the tube in meters
layout = 't' # t/s triangle or square
k_wall = 15  # thermal conductivity of tube wall in W/mK
# shell side parameters (ketone)
m_dos = 22.22
c_pShell = 2200.8 #specific heat capacity in J/kgK (ketone)
muShell = 0.4131
tempInShell = 283.15
tempOutShell = 303.15
sgShell = 810 # density of the shell's fluid in kg/m³
D_shell = 0.5  # shell diameter in meters
baffle_cut = 0.25  # baffle cut as a fraction of the shell diameter
Q = 50000  # total heat transfer in Watts
deltaT_lm = tlm(tempInShell, tempOutShell, tempInTube, tempOutTube)  # log mean temperature difference in K

# Calculate coefficients
h_tube = tube_side_coefficient(d_i, m_dot, c_p, mu, k, l)
h_shell = shell_side_coefficient(d_o, baffle_cut, d_o*1.25, m_dos, c_pShell, muShell, k, l, sgShell, D_shell)
U = overall_heat_transfer_coefficient(h_tube, h_shell, d_o, d_i, k_wall)
A = required_area(Q, U, deltaT_lm)
nTube = numOfTubes(D_shell, d_o, layout)
#length, nTube = estimate_tube_number(A, d_o, l, nTube)

print(f"Tube-side heat transfer coefficient: {h_tube:.2f} W/m²K")
print(f"Shell-side heat transfer coefficient: {h_shell:.2f} W/m²K")
print(f"Overall heat transfer coefficient: {U:.2f} W/m²K")
print(f"Required heat transfer area: {A:.2f} m²")
