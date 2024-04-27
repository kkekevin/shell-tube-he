import numpy as np

def calculate_pressure_drop(length, diameter, flow_rate, density, viscosity):
    # Constants
    g = 9.81  # acceleration due to gravity in m/s^2
    
    # Calculate velocity
    area = np.pi * (diameter / 2)**2
    velocity = flow_rate / area

    # Calculate Reynolds number
    Re = density * velocity * diameter / viscosity
    
    # Calculate friction factor using Churchill's equation (applicable over all regimes)
    A = (2.457 * np.log(1 / ((7 / Re)**0.9 + 0.27 * (0.0001/diameter))))**16
    B = (37530 / Re)**16
    f = 8 * ((8 / Re)**12 + 1 / (A + B)**1.5)**(1/12)
    
    # Calculate pressure drop using Darcy-Weisbach equation
    pressure_drop = f * (length / diameter) * (density * velocity**2 / 2) / g
    
    return pressure_drop

# Example parameters (SI Units)
length = 10.0       # length of the tube in meters
diameter = 0.025    # diameter of the tube in meters
flow_rate = 0.1     # flow rate in m^3/s
density = 1000      # density of water in kg/m^3
viscosity = 0.001   # viscosity of water in Pa.s

pressure_drop = calculate_pressure_drop(length, diameter, flow_rate, density, viscosity)
print(f"Pressure Drop: {pressure_drop:.2f} Pa")
