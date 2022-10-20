The Atmospheric Radiative Transfer Simululator [ARTS](https://radiaitvetransfer.org) is one of several widely-available
line-by-line models. ARTS is forward-looking in having developed a Python interface we can exploit in these exercises. 

### Installing the line-by-line code and absorption data 

The `pyarts_utilities.py` code in this directory wraps two modules adapted from [lectures](https://github.com/atmtools/arts-lectures) 
given at the University of Hamburg by Manfred Brath, Oliver Lemke, Stefan Buehler, and their colleagues. The lecture 
modules in turn wrap calls to `pyarts`. You'll need to install `pyarts` and its dependencies including `typhon` (utilities
developed by the same group). 

I suggest using a conda-based package manager -- I like [mamba](https://mamba.readthedocs.io/) -- to create a virtual environment. 
`pyarts` and `typhon` are both available through the `rttools` channel, or you may use the `environment.yml` file in this directory, i.e. 
```
mamba env create -f environment.yml
``` 
You'll likely want to add other plotting and analysis software to this environment to make plots, etc. 

You will also need to download the current snapshot of the absorption [catalog data](https://www.radiativetransfer.org/misc/download/unstable/) (the 'cat' file, not the 'xml' file). You will need to provide the downloaded location to the codes. 

### An extremely simple pyARTSinterface for two tasks

ARTS is general and flexible and can be somewhat involved to use. Module `pyarts_utils` includes three functions: 
- `calculate_absxsec_wn()` computes the spectrally-resolved absorption cross-section ($m^2$ per molecule) across a 
user-specfied range of wavenumbers for a user-specified set of gases.
- `calc_olr_wn()` computes the spectrally-resolved top-of-atmosphere outgoing longwave radiation. 
- `create_arts_atm()` creates the representation of the atmosphere required by `calc_olr_wn()` from vectors 
of temperature and pressure and a dictionary of gas volume mixing rates (keys are chemical formulae like `"H2O"`, values are 
volume mixing ratios)

Both computational functions return wavenumber and the desired result (i.e. 'wn, beta = calculate_absxsec_wn()'). Both need the 
root location of the absorption catalog supplied as a `pathlib.Path()` in argument `basefile.` The default spectral resolution is relatively 
coarse. 

The module also include variables that can be used to specify the absorption features used in computations: 
- variables `H2O,CO2,CH4,N2O,O3` specify line absorption by each gas
- `H2O_self, H2O_foreign` specify the self- and foreign continuum for water vapor (MT_CKD 3.5 in case this matters)
- `H2O_plus, CO2_plus` include both lines and coninuua. 