# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the GNU License

import click

"""
Command line interface for the SL3ME package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.1"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Oct 2018"

# This is used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    CLI tools for the SL3ME package.
    """
    pass

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("data_file", nargs=1)
@click.option("--direct_band_gap", "-D", default=0.0)
@click.option("--indirect_band_gap", "-I", default=0.0)
@click.option("--thickness", "-t", default=5e-7)
@click.option("--temperature", "-T", default=300)
def calculate(data_file, direct_band_gap, indirect_band_gap, thickness,
              temperature):
    """
    Calculate the SLME from the optical data.
    """
    from SL3ME.cli.commands import calculate_slme

    calculate_slme(data_file=data_file,
                   direct_band_gap=direct_band_gap,
                   indirect_band_gap=indirect_band_gap,
                   thickness=thickness,
                   temperature=temperature)