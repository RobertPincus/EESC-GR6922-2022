from pyarts_utils      import calculate_absxsec_wn, H2O_plus, calc_olr_wn
from pathlib import Path 
import numpy as np

wn, beta = calculate_absxsec_wn(H2O_plus, 
                                wn_min = 800, wn_max = 1000, wn_num = 1000, 
                                arts_data_root = Path("/Users/robert/Codes/arts-cat-data"))

wn, olr  = calc_olr_wn(arts_data_root = Path("/Users/robert/Codes/arts-cat-data"))
print(np.trapz(olr, wn))