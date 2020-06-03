##  En este script se calculan los acumulados de precipitacion
##  para cada metodo.
##  Como entrada toma los datos de las tasas de precipitacion
##  que fueron obtenidas previamente.

import pyart
from numpy import isnan
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib as mpl
from numpy import array
from os import listdir

dir_images='/home/romina/radar/doctorado/Datos/20100118/datos_PP_00654/'

fecha='20100118'

#variable='rain_kdp_gu'
variable='rain_att'
#variable='rain_dbz'

var='_att_'

files = sorted([f for f in listdir(dir_images) if f.endswith('rain_att.nc') ])

the_maximums=0.

for fi in files:
    img = pyart.io.read(dir_images + fi)
    
    print(fi)
    
    dato=img.fields[variable]['data']
    
    dato[isnan(dato)]=0.
    dato[dato<0.]=0.
    
    
    desaz=0
    sweep=0

    i1_1 = img.sweep_start_ray_index['data'][sweep]
    i2_1 = img.sweep_end_ray_index['data'][sweep]
    azimuth_1=img.azimuth['data'][i1_1:i2_1]
    azdiff_1=abs(azimuth_1-desaz)
    ray_num_1 = azdiff_1.argmin()+i1_1
    
    aux_a = dato[:ray_num_1]
    aux_b = dato[ray_num_1:]
    dato_flip = array(list(aux_b)+list(aux_a))
    
    the_maximums = the_maximums + dato_flip
    
acumulado=the_maximums/6.

img.add_field_like(variable,'acumulado',acumulado,replace_existing='True')    
img.azimuth['data'] = np.array(range(0,361))

pyart.io.write_cfradial(dir_images + fecha + var + 'acumulado.nc',img) 
