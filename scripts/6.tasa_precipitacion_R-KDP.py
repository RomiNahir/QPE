# Este script realiza estimaciones de precipitaciones con kdp

import numpy as np
import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import listdir
from sys import argv
from funciones import VariableUnica


HOME='/home/romina/radar/doctorado/'

# el path hay que indicar la facha que se quiere utilizar

path=argv[1]

for fecha in list(set([f[0:12] for f in listdir(HOME + path+'/datos_PHIDP') if f.endswith('.nc')])):
    print(fecha)
    anio=fecha[0:4]
    mes=fecha[4:6]
    dia=fecha[6:8]
    hora=fecha[8:12]

    RADAR=pyart.io.read(HOME + path + 'datos_PHIDP/' + fecha + 'cor.nc')
    RADARkdp=pyart.io.read(HOME + path + 'datos_KDP/' + fecha + 'kdp.nc')
   
    zdr = RADAR.fields['ZDR']['data']
    dbz = RADAR.fields['dBZ']['data']
    rho = RADAR.fields['RhoHV']['data']
    phi = RADAR.fields['maskPhi']['data']
    kdp = RADARkdp.fields['mask_kdp']['data']
    
    nb=RADAR.ngates
    nr=RADAR.nrays
    
    ######### Precipitacion con kdp Gu et al (2011) ##############

    a_kdp_gu=25.1
    b_kdp_gu=0.777

    rain_kdp_gu=np.zeros((nr,nb))
    
    rain_kdp_gu = a_kdp_gu*np.ma.power(kdp, b_kdp_gu)

    #rain_kdp_gu=a_kdp_gu*kdp**b_kdp_gu
        
    mask_rain_kdp_gu=np.ma.masked_where(np.isnan(RADAR.fields['mask_ref']['data']),rain_kdp_gu)

    aa_kdp_gu=np.ma.filled(mask_rain_kdp_gu,fill_value=np.nan)
    bb_kdp_gu = np.ma.masked_invalid(aa_kdp_gu)
    
#### Agrego las variables nuevas ####

    RADAR.add_field_like('uncorrected_differential_phase','pp_kdp',bb_kdp_gu)
    
    RAIN=VariableUnica(RADAR, bb_kdp_gu, 'dBZ','rain_kdp_gu')
    
    pyart.io.write_cfradial(HOME + path + 'datos_PP_KDP/' + fecha+'rain_kdp.nc',RAIN)
