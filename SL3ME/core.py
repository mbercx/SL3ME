# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

import numpy as np
import matplotlib as plt

import scipy.constants as constants

from scipy.integrate import simps
from scipy.optimize import minimize

"""
Core classes and methods of the SL3ME package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "May 2018"


def slme(material_energy_for_absorbance_data,
         material_absorbance_data,
         material_direct_allowed_gap,
         material_indirect_gap, thickness=50E-6, temperature=293.15,
         absorbance_in_inverse_centimeters=False,
         cut_off_absorbance_below_direct_allowed_gap=True,
         plot_current_voltage=False):
    """
    Calculate the

    IMPORTANT NOTES:
    1) Material calculated absorbance is assumed to be in m^-1, not cm^-1!
        (Most sources will provide absorbance in cm^-1, so be careful.)

    2) The default is to remove absorbance below the direct allowed gap.
        This is for dealing with broadening applied in DFT absorbance
        calculations. Probably not desired for experimental data.

    3) We can calculate at different temperatures if we want to, but 25 C /
        293.15 K is the standard temperature assumed if not specified

    4) If absorbance is in cm^-1, multiply values by 100 to match units
        assumed in code

    Args:
        material_energy_for_absorbance_data:
        material_absorbance_data:
        material_direct_allowed_gap:
        material_indirect_gap:
        thickness:
        temperature:
        absorbance_in_inverse_centimeters:
        cut_off_absorbance_below_direct_allowed_gap:
        plot_current_voltage:

    Returns:
        The calculated maximum efficiency.

    """

    # Defining constants for tidy equations
    c = constants.c  # speed of light, m/s
    h = constants.h  # Planck's constant J*s (W)
    h_e = constants.h / constants.e  # Planck's constant eV*s
    k = constants.k  # Boltzmann's constant J/K
    k_e = constants.k / constants.e  # Boltzmann's constant eV/K
    e = constants.e  # Coulomb

    # Make sure the absorption coefficient has the right units (m^{-1})
    if absorbance_in_inverse_centimeters:
        material_absorbance_data = material_absorbance_data * 100

    # Load the Air Mass 1.5 Global tilt solar spectrum
    solar_spectrum_data_file = os.path.join(
        os.path.join(os.path.dirname(__file__), os.pardir),
        "data",
        "am1.5G.dat"
    )

    solar_spectra_wavelength, solar_spectra_irradiance = np.loadtxt(
        solar_spectrum_data_file, usecols=[0, 1], unpack=True, skiprows=2
    )

    solar_spectra_wavelength_meters = solar_spectra_wavelength * 1e-9

    delta = material_direct_allowed_gap - material_indirect_gap
    fr = np.exp(-delta / (k_e * temperature))

    # need to convert solar irradiance from Power/m**2(nm) into
    # photon#/s*m**2(nm) power is Watt, which is Joule / s
    # E = hc/wavelength
    # at each wavelength, Power * (wavelength(m)/(h(Js)*c(m/s))) = ph#/s
    solar_spectra_photon_flux = solar_spectra_irradiance * (
        solar_spectra_wavelength_meters / (h * c))

    ### Calculation of total solar power incoming
    power_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    ### calculation of blackbody irradiance spectra
    ## units of W/(m**3), different than solar_spectra_irradiance!!! (This
    # is intentional, it is for convenience)
    blackbody_irradiance = (2.0 * h * c ** 2 /
                            (solar_spectra_wavelength_meters ** 5)) \
                           * (1.0 / ((np.exp(h * c / (
        solar_spectra_wavelength_meters * k * temperature))) - 1.0))

    # I've removed a pi in the equation above - Marnik Bercx

    ## now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (
        solar_spectra_wavelength_meters / (h * c))

    ## units of nm
    material_wavelength_for_absorbance_data = ((c * h_e) / (
        material_energy_for_absorbance_data + 0.00000001)) * 10 ** 9

    ### absorbance interpolation onto each solar spectrum wavelength
    from scipy.interpolate import interp1d
    ## creates cubic spline interpolating function, set up to use end values
    #  as the guesses if leaving the region where data exists
    material_absorbance_data_function = interp1d(
        material_wavelength_for_absorbance_data, material_absorbance_data,
        kind='cubic',
        fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
        bounds_error=False
    )

    material_interpolated_absorbance = np.zeros(
        len(solar_spectra_wavelength_meters))
    for i in range(0, len(solar_spectra_wavelength_meters)):
        ## Cutting off absorption data below the gap. This is done to deal
        # with VASPs broadening of the calculated absorption data


        if solar_spectra_wavelength[i] < 1e9 * ((c * h_e) /
                                                    material_direct_allowed_gap) \
                or cut_off_absorbance_below_direct_allowed_gap == False:
            material_interpolated_absorbance[
                i] = material_absorbance_data_function(
                solar_spectra_wavelength[i])

    absorbed_by_wavelength = 1.0 - np.exp(-2.0 *
                                          material_interpolated_absorbance
                                          * thickness)
    # np.ones(material_interpolated_absorbance.size)
    print(str(type(absorbed_by_wavelength)))

    ###  Numerically integrating irradiance over wavelength array
    ## Note: elementary charge, not math e!  ## units of A/m**2   W/(V*m**2)
    J_0_r = e * np.pi * simps(blackbody_photon_flux * absorbed_by_wavelength,
                              solar_spectra_wavelength_meters)

    J_0 = J_0_r / fr

    ###  Numerically integrating irradiance over wavelength array
    ## elementary charge, not math e!  ### units of A/m**2   W/(V*m**2)
    J_sc = e * simps(solar_spectra_photon_flux * absorbed_by_wavelength,
                     solar_spectra_wavelength)

    #    J[i] = J_sc - J_0*(1 - exp( e*V[i]/(k*T) ) )
    #   #This formula from the original paper has a typo!!
    #    J[i] = J_sc - J_0*(exp( e*V[i]/(k*T) ) - 1)
    #   #Bercx chapter and papers have the correct formula (well,
    #   the correction on one paper)
    def J(V):
        J = J_sc - J_0 * (np.exp(e * V / (k * temperature)) - 1.0)
        return J

    def power(V):
        p = J(V) * V
        return p

    def neg_power(V):  # Minimizing the negative of power to optimize power
        return -power(V)

    results = minimize(neg_power, x0=np.array([0.0000001]))
    # passing a function as a variable, preconditioning at
    # V=0, since this is a smooth function results.x are the inputs as an
    # array that give the highest value of Power

    # I've changed the preconditioning to a very small value, as per Kamal
    # Choudhary's suggestion - Marnik Bercx
    # TODO figure out how that minimize method works

    V_Pmax = float(results.x)
    P_m = power(V_Pmax)

    efficiency = P_m / power_in

    # This segment isn't needed for functionality at all, but can display a
    # plot showing how the maximization of power by choosing the optimal
    # voltage value works
    if plot_current_voltage:
        V = np.linspace(0, 2, 200)
        plt.plot(V, J(V))
        plt.plot(V, power(V), linestyle='--')
        plt.show()
        print(power(V_Pmax))

    return efficiency
