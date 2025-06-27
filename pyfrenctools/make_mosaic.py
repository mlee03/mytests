import ctypes
import numpy as np
import numpy.typing as npt
from typing import Union, List
import xarray as xr

from make_mosaic_utils import (
    AtmComponent, LndComponent, OcnComponent, set_component
)
from mosaicobj import MosaicObj

#atm 
atm_mosaic = MosaicObj(input_dir="./input", mosaic_file="C48_mosaic.nc").read()
atm_mosaic.get_grid()

nxy = atm_mosaic.grid['tile1'].nx * atm_mosaic.grid['tile1'].ny
atm_masks = [None] * atm_mosaic.ntiles
for imask in atm_masks: imask = np.ones(nxy, dtype=np.float64)    

atm, atm_argtype = set_component(mosaic=atm_mosaic, Component=AtmComponent, mask=atm_masks)


#lnd
lnd_mosaic = MosaicObj(input_dir="./input", mosaic_file="C48_mosaic.nc").read()
lnd_mosaic.get_grid()
nxy = lnd_mosaic.grid['tile'].nx * lnd_mosaic.grid['tile'].ny
lnd_masks = [None] * lnd_mosaic.ntiles
for imask in lnd_masks: imask = np.ones(nxy, dtype=np.float64)

lnd, lnd_argtype = set_component(mosaic=lnd_mosaic, Component=LndComponent, mask=lnd_masks)

#ocn
ocn_mosaic = MosaicObj(input_dir="./input", mosaic_file="ocean_mosaic.nc").read()
ocn_mosaic.get_grid()

ocn, ocn_argtype = set_component(mosaic=ocn_mosaic, Component=OcnComponent)                    

exit()

lib = ctypes.cdll.LoadLibrary("./c_install/clib.so")
lib.make_coupler_mosaic.restype = ctypes.c_int
lib.make_coupler_mosaic.argtypes = [ctypes.c_int, #ntile_lnd
                                    ctypes.c_int, #ntile_lnd
                                    ctypes.c_int, #ntile_ocn
                                    ctypes.c_int, #ocn_south_ext
                                    ctypes.POINTER(atm_argtype), #lnd
                                    ctypes.POINTER(lnd_argtype), #lnd
                                    ctypes.POINTER(ocn_argtype)  #ocn
]
    
lib.make_coupler_mosaic(atm_mosaic.ntiles,
                        lnd_mosaic.ntiles,
                        ocn_mosaic.ntiles,
                        1, ctypes.byref(atm), ctypes.byref(lnd), ctypes.byref(ocn))
