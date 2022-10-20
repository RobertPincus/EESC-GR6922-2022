""" Utilities for using pyARTS 
"""
from pathlib import Path 
from scipy.constants   import speed_of_light as c
from absorption_module import calculate_absxsec
from olr_module        import calc_olr, Change_T_with_RH_const
from pyarts            import xml

H2O,CO2,CH4,N2O,O3 = 'H2O','CO2','CH4','N2O','O3'
H2O_self = "H2O-SelfContCKDMT350"
H2O_foreign = "H2O-ForeignContCKDMT350"
H2O_plus = 'H2O, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350'
CO2_plus = 'CO2, CO2-CKDMT252'
##################################################################
def spectral_grid_from_wavenumbers(wn_min=0, wn_max=3200, 
	                               wn_spacing=.05, wn_num=None): 
	assert wn_max > wn_min
	fmin = wn_min * c * 100
	fmax = wn_max * c * 100
	fnum = int((wn_max - wn_min)/wn_spacing + 1)
	if  wn_num is not None: fnum = wn_num
	return ( {"fmin":fmin, "fmax":fmax, "fnum":fnum} )
##################################################################
def calculate_absxsec_wn(species=H2O_plus,
                      pressure=800e2,
                      temperature=300.0,
                      wn_min=10e9,
                      wn_max=2000e9,
                      wn_num=10001,
                      vmr=0.05,
                      arts_data_root=Path("./arts-cat-data")):
	""" Wrapper for absorption_module.py:calculate_absxsec from ARTS lectures 
    Parameters:
        species (str): Absorption species name.
        pressure (float): Atmospheric pressure [Pa].
        temperature (float): Atmospheric temperature [K].
        wn_min (float): Minimum wavenumber [cm-1].
        wn_max (float): Maximum frequency [cm-1].
        wn_num (int): Number of frequency grid points.
        vmr (float): Volume mixing ratio. This is mainly important for the
                     water vapor continua.
        arts_data_root (pathlib.Path): Root location of arts-cat-data

    Returns:
        ndarray, ndarray: Wavenumber grid [cm-1], Abs. cross sections [m^2]
	"""
	fgrid = spectral_grid_from_wavenumbers(wn_min, wn_max, wn_num=wn_num)
	freq, beta = calculate_absxsec(species=species,
					  pressure=pressure,
                      temperature=temperature,
                      vmr=vmr,
                      basename = str(arts_data_root.joinpath("lines")) + "/", 
                      **fgrid)
	return (freq/(c * 100), beta)

##################################################################
def calc_olr_wn(wn_min=1e-9, wn_max=3250., wn_num=10001,
	            delta_T=0.,  CO2_scale=1., 
	            arts_data_root=Path("./arts-cat-data")): 
	""" Wrapper for olr_module.py:calc_olr from ARTS lectures 
	Calculates spectrally-resolved OLR for the (potentially-modified) AFGL mid-latitude summer atmosphere

    Parameters:
        wn_min (float): Minimum wavenumber [cm-1].
        wn_max (float): Maximum frequency [cm-1].
        wn_num (int): Number of frequency grid points.
        delta_T (float): temperature change uniformly added to profiles; relative humidity is conserved
        CO2_scale (float): multiplier for CO2 concentrations
        arts_data_root (pathlib.Path): Root location of arts-cat-data

    Returns:
        ndarray, ndarray: Wavenumber grid [cm-1], OLR[W/m^2-cm^-1]
    """
	assert wn_min >= 1e-9
	atmfield = xml.load("input/midlatitude-summer.xml")
	atmfield.set("T", atmfield.get("T") + 0)


	# Scale the CO2 concentration
	if CO2_scale is not None: 
		atmfield.scale("abs_species-CO2", CO2_scale)

	# Add a constant value to the temperature  
	# with constant  relative humidity 
	if delta_T is not None: 
		atmfield = Change_T_with_RH_const(atmfield, DeltaT=delta_T)

	fgrid = spectral_grid_from_wavenumbers(wn_min, wn_max, wn_num=wn_num)
	freq, olr = calc_olr(atmfield,
				         basename = str(arts_data_root.joinpath("lines")) + "/", 
				         **fgrid)
	return (freq/(c * 100), olr * (c * 100)) 

##################################################################
