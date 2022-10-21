from pyarts_utils import calculate_absxsec_wn, H2O_plus, create_arts_atm, calc_olr_wn
from pathlib import Path
import numpy as np

wn, beta = calculate_absxsec_wn(
    H2O_plus,
    wn_min=800,
    wn_max=1000,
    wn_num=1000,
    arts_data_root=Path("/Users/robert/Codes/arts-cat-data"),
)


#
# Absurd profiles of T, z, q, and CO2
#
pressure = np.logspace(np.log10(101300), np.log10(10000), 20)
T = np.linspace(290, 200, len(pressure))
vmr_h2o = 0.01 * np.exp(-np.linspace(0, 16000, len(pressure)) / 2000)
vmr_co2 = np.ones(np.size(pressure)) * 400 / 1e6
vmr_co2 = vmr_co2[0]
atm = create_arts_atm(pressure, T, {"h2o": vmr_h2o, "co2": vmr_co2})

wn, olr = calc_olr_wn(atm, arts_data_root=Path("/Users/robert/Codes/arts-cat-data"))
print("Integrated OLR is ", np.trapz(olr, wn))
