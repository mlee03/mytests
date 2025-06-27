import ctypes
import numpy as np
import numpy.typing as npt
from typing import Union, List
import xarray as xr

from mosaicobj import MosaicObj

class AtmComponent(ctypes.Structure): pass
class LndComponent(ctypes.Structure): pass
class OcnComponent(ctypes.Structure): pass

def set_component(mosaic: MosaicObj, Component: Union[AtmComponent, LndComponent, OcnComponent],
                  mask: List[npt.NDArray[np.float64]] = None, area: List[npt.NDArray[np.float64]] = None):    
    
    Component._fields_ = [("itile", ctypes.c_int),
                          ("nx",  ctypes.c_int),
                          ("ny", ctypes.c_int),
                          ("x", ctypes.POINTER(ctypes.c_double)),
                          ("y", ctypes.POINTER(ctypes.c_double)),
                          ("mask", ctypes.POINTER(ctypes.c_double)),
                          ("area", ctypes.POINTER(ctypes.c_double))
    ]
    
    componentss = Component * mosaic.ntiles
    components = componentss()

    for itile in range(mosaic.ntiles):
        tile = mosaic.gridtiles[itile]
        igrid = mosaic.grid[tile]
        components[itile].itile = itile
        components[itile].nx = igrid.nx
        components[itile].ny = igrid.ny
        components[itile].x = igrid.x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        components[itile].y = igrid.y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        components[itile].mask = mask if mask is None else mask[itile].data_as(ctypes.POINTER(ctypes.c_double))
        components[itile].area = area if area is None else area[itile].data_as(ctypes.POINTER(ctypes.c_double))
    return components, componentss


def extend_south(ocn_mosaic: MosaicObj):

    tiny_value = 1.e-7
    min_atm_lat = np.radians(-90.0)
    
    if ocn_mosaic.grids['tile'].y[0][0] > min_atm_lat + tiny_value:
        #extend
        x = np.



    
