# Este script realiza estimaciones de precipitaciones con
# reflectividad corregida


import numpy as np
import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import listdir
from sys import argv
from funciones import VariableUnica

HOME='/home/romina/radar/doctorado/'

# el path hay que indicar la fecha que se quiere utilizar

path=argv[1]

for fecha in list(set([f[0:12] for f in listdir(HOME + path+'/datos_PHIDP') if f.endswith('cor.nc')])):
    print(fecha)
    
    anio=fecha[0:4]
    mes=fecha[4:6]
    dia=fecha[6:8]
    hora=fecha[8:12]
    
    RADAR=pyart.io.read(HOME + path + 'datos_ATT_00654/' + fecha + 'att.nc')
    RADARatt=pyart.io.read(HOME + path + 'datos_ATT_00654/' + fecha + 'att.nc')
        
    dbz = RADAR.fields['mask_zcor']['data']
    
    nb=RADAR.ngates
    nr=RADAR.nrays
    
   ##### Precipitacion con Z #######
    
    a_dbz=300
    b_dbz=1.4
    
    ## Transformo el dbz en mm6 m-3

    dbz_lineal = 10.0 ** (0.1 * dbz)
    rain_dbz=np.zeros((nr,nb))

    rain_dbz=10**(np.log10(dbz_lineal/a_dbz)/b_dbz)
    
    #rain_dbz=10*np.ma.power((dbz_lineal/a_dbz,1))/b_dbz   
    
    #mask_rain_dbz=np.ma.masked_where(np.isnan(RADAR.fields['mask_ref']['data']),rain_dbz)
    #aa_dbz_2=np.ma.filled(mask_rain_dbz,fill_value=np.nan)
    #bb_dbz = np.ma.masked_invalid(aa_dbz)

    mask_rain_dbz=np.ma.masked_where(np.isnan(RADARatt.fields['mask_att']['data']),rain_dbz)
    aa_dbz=np.ma.filled(mask_rain_dbz,fill_value=np.nan)
    bb_dbz = np.ma.masked_invalid(aa_dbz)


#### Agrego las variables nuevas ####

    RADARatt.add_field_like('mask_att','pp_dbz',bb_dbz)

    RAIN=VariableUnica(RADARatt, bb_dbz, 'mask_att','rain_dbz')
    
    pyart.io.write_cfradial(HOME + path +'datos_PP_ZR/' +fecha+'rain_dbz_cor.nc',RAIN)
