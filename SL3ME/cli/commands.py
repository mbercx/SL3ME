# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the GNU License

import numpy as np

from SL3ME.core import slme

"""
Scripts for the CLI commands.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Oct 2018"


def calculate_slme(data_file, direct_band_gap, indirect_band_gap, thickness,
                   temperature):
    """
    Calculate the SLME from an optical data file.

    Args:
        data_file:
        direct_band_gap:
        indirect_band_gap:
        thickness:
        temperature:

    Returns:

    """
    # Load the optical data
    energy, absorbance = np.loadtxt(data_file, usecols=[0, 1], unpack=True)


    efficiency = slme(
        material_energy_for_absorbance_data=energy,
        material_absorbance_data=absorbance,
        material_direct_allowed_gap=direct_band_gap,
        material_indirect_gap=indirect_band_gap,
        thickness=thickness,
        temperature=temperature
    )

    print("The calculated efficiency = " + str(efficiency*100) + " %")