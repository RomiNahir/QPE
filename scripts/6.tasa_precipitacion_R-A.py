# Este script realiza estimaciones de precipitaciones con atenuacion

import numpy as np
import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import listdir
from sys import argv
from funciones import VariableUnica

HOME='/home/romina/radar/doctorado/'

# el path hay que indicar la facha que se quiere utilizar /Datos/fecha/

path=argv[1]

for fecha in list(set([f[0:12] for f in listdir(HOME + path+'/datos_PHIDP') if f.endswith('.nc')])):
    print(fecha)
    anio=fecha[0:4]
    mes=fecha[4:6]
    dia=fecha[6:8]
    hora=fecha[8:12]

    RADAR=pyart.io.read(HOME + path + 'datos_PHIDP/' + fecha + 'cor.nc')
    RADARkdp=pyart.io.read(HOME + path + 'datos_KDP/' + fecha + 'kdp.nc')
    RADARatt=pyart.io.read(HOME + path + 'datos_ATT_008/metodo_simple/' + fecha + 'att.nc')
    
    att = RADARatt.fields['mask_att']['data']

    att[att>10]=np.nan
    att = np.ma.masked_invalid(att)    

    nb=RADAR.ngates
    nr=RADAR.nrays
    
    a_att=294.0
    b_att=0.89
    
    rain_att=np.zeros((nr,nb))
 
    rain_att = a_att*np.ma.power(att, b_att)

            
    mask_rain_att=np.ma.masked_where(np.isnan(RADAR.fields['mask_ref']['data']),rain_att)
    
    aa_att=np.ma.filled(mask_rain_att,fill_value=np.nan)
    bb_att=np.ma.masked_invalid(aa_att)
    

#### Agrego las variables nuevas ####

    RADAR.add_field_like('uncorrected_differential_phase','pp_att',bb_att)

    
    RAIN=VariableUnica(RADAR, bb_att, 'dBZ','rain_att')
    
    pyart.io.write_cfradial(HOME + path +'datos_PP_008/metodo_simple/' + fecha + 'rain_att.nc',RAIN)
