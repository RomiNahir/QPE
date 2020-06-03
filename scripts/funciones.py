# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 11:13:07 2016

@author: romina
"""
import numpy as np
from numpy import convolve
import copy
import pyart
 
def promediomovil(values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

def VariableUnica(objRADAR, datos, var_copiar, var_nueva):
    """Crea un nuevo objeto Radar haciendo una copia del que le pasamos 
    y agrega una nueva variable similar a alguna que exista y borra el resto
    Parámetros:
        - objRADAR: 
        -datos:
        ... """
    copia = copy.deepcopy(objRADAR)
    copia.add_field_like(var_copiar, var_nueva, datos)
    
    # Borra el resto de las variables
    #for k in copia.fields.keys():
    for k in list(copia.fields):
        if k != var_nueva:
            copia.fields.pop(k)

    return copia

def AlturaHaz(radar,sweep):

    ranges = radar.range['data']
    elevation = radar.fixed_angle['data'][sweep]
    radar_height = radar.altitude['data']
    Re = 6371.0 * 1000.0
    p_r = 4.0 * Re / 3.0

    z = radar_height + (ranges ** 2 + p_r ** 2 + 2.0 * ranges * p_r *
                        np.sin(elevation * np.pi / 180.0)) ** 0.5 - p_r
    
    return z                    
    
    
def getCartesianGrid(radar, data):

        metros = 1000.0
                  
        stopRange = 240
        z = 1
        grid = pyart.map.grid_from_radars(
            (radar,),
            grid_shape=(z, 480, 480),
            grid_limits=((0, z ), 
                         (-stopRange * metros, stopRange * metros),
                         (-stopRange * metros, stopRange * metros)),
            weighting_function='Barnes2',
            fields=[data],
            min_radius=1100,
            grid_origin=(radar.latitude['data'][0], radar.longitude['data'][0])

        )

        return grid
        
def saveToGTiff(grid, field, outFilePath, outFileName):
    pyart.io.write_grid_geotiff(
        grid,
        filename=outFilePath + outFileName + ".tif",
        field=field,
        #level=0,
        #rgb=True,
        #cmap= self.__rainbowRadar.getRadarVariable()[2],
        vmin= 0.,#self.__rainbowRadar.getRadarVariable()[3],
        vmax= 1.,#self.__rainbowRadar.getRadarVariable()[4],
        warp=True
        )            
