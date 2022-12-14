"""Calculate and plot absorption cross sections.
   Taken from a snapshot of https://github.com/atmtools/arts-lectures compatible 
   with pyARTS 2.6.4. Codes by Manfred Brath, Oliver Lemke, and colleagues 
"""
import re

import numpy as np
import pyarts
from pathlib import Path


def tag2tex(tag):
    """Replace all numbers in a species tag with LaTeX subscripts."""
    return re.sub("([a-zA-Z]+)([0-9]+)", r"\1$_{\2}$", tag)


def calculate_absxsec(
    species="H2O, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350",
    pressure=800e2,
    temperature=300.0,
    fmin=10e9,
    fmax=2000e9,
    fnum=10_000,
    lineshape="LP",
    normalization="RQ",
    verbosity=0,
    vmr=0.05,
    basename="lines/",
    lines_off=0,
):
    """Calculate absorption cross sections.

    Parameters:
        species (str): Absorption species name.
        pressure (float): Atmospheric pressure [Pa].
        temperature (float): Atmospheric temperature [K].
        fmin (float): Minimum frequency [Hz].
        fmax (float): Maximum frequency [Hz].
        fnum (int): Number of frequency grid points.
        lineshape (str): Line shape model.
                            Available options:
                            DP        -      Doppler profile,
                            LP        -      Lorentz profile,
                            VP        -      Voigt profile,
                            SDVP      -      Speed-dependent Voigt profile,
                            HTP       -      Hartman-Tran profile.
        normalization (str): Line shape normalization factor.
                            Available options:
                            VVH       -      Van Vleck and Huber,
                            VVW       -      Van Vleck and Weisskopf,
                            RQ        -      Rosenkranz quadratic,
                            None      -      No extra normalization.
        verbosity (int): Set ARTS verbosity (``0`` prevents all output).
        vmr (float): Volume mixing ratio. This is mainly important for the
                     water vapor continua.
        arts_data_root (pathlib.Path): Root location of arts-cat-data
        lines_off (int): Switch off lines, if no contnua is selected,
                         absorption will be zero.

    Returns:
        ndarray, ndarray: Frequency grid [Hz], Abs. cross sections [m^2]
    """
    # Create ARTS workspace and load default settings
    ws = pyarts.workspace.Workspace(verbosity=0)
    ws.LegacyContinuaInit()
    ws.water_p_eq_agendaSet()
    ws.PlanetSet(option="Earth")

    ws.verbositySetScreen(ws.verbosity, verbosity)

    # We do not want to calculate the Jacobian Matrix
    ws.jacobianOff()

    # Define absorption species
    ws.abs_speciesSet(species=[species])
    ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename=basename)
    ws.abs_lines_per_speciesLineShapeType(option=lineshape)
    ws.abs_lines_per_speciesCutoff(option="ByLine", value=750e9)
    ws.abs_lines_per_speciesNormalization(option=normalization)
    if lines_off:
        ws.abs_lines_per_speciesSetEmpty()

    # Create a frequency grid
    ws.VectorNLinSpace(ws.f_grid, fnum, fmin, fmax)

    # Throw away lines outside f_grid
    ws.abs_lines_per_speciesCompact()

    # Atmospheric settings
    ws.AtmosphereSet1D()
    ws.stokes_dim = 1

    # Setting the pressure, temperature and vmr
    ws.rtp_pressure = float(pressure)  # [Pa]
    ws.rtp_temperature = float(temperature)  # [K]
    ws.rtp_vmr = np.array([vmr])  # [VMR]
    ws.Touch(ws.rtp_nlte)

    # isotop
    ws.isotopologue_ratiosInitFromBuiltin()

    # Calculate absorption cross sections
    ws.lbl_checkedCalc()
    ws.propmat_clearsky_agenda_checked = 1
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddLines()
    ws.propmat_clearskyAddConts()

    # Convert abs coeff to cross sections on return
    number_density = pressure * vmr / (pyarts.arts.constant.k * temperature)

    return (
        ws.f_grid.value.value.copy(),
        ws.propmat_clearsky.value.data.value[0, 0, :, 0].copy() / number_density,
    )
