# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:52:45 2017

@author: romina
"""



### Este script calcula la atenuacion con el metodo simple
### de A=alpha*kdp**b y la correcion en reflectividad
### Es importante que previamente sepa cual es el alpha
### que se va a usar porque puede cambiar según la situación

#### Esta version pone como mascara al kdp porque asi se van
#### pixels con datos erroneos en el centro que son aislados
#### esto hay que corregirlo en un futuro con algun filtro textural

import numpy as np
import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
import rango as ran
from scipy.integrate import cumtrapz
from sistema import smooth_and_trim
import numpy.ma as ma
import glob, os
from sys import argv
from funciones import VariableUnica

HOME='/home/romina/radar/doctorado/Datos/'


# el path es /fecha/datos_PHIDP/ donde estan los corregidos de phidp

path=argv[1]


os.chdir(HOME+path)
for file in glob.glob("*cor.nc"):
    print(file)
    
    RADAR = pyart.io.read(HOME + path + file)
    
    RADARkdp=pyart.io.read(HOME + '20090721/' + 'datos_KDP/' + file[0:12] + 'kdp.nc')

    
    levels=np.unique(RADAR.elevation['data'])

    #los coeficientes son para banda C calculados con el alpha_beta

    b=0.8    #coeficiente b
    alpha=0.08  #coeficiente alfa
    
    rho = RADAR.fields['RhoHV']['data']
    maskphi=RADAR.fields['maskPhi']['data']
    maskref=RADAR.fields['mask_ref']['data']
    kdp=RADARkdp.fields['mask_kdp']['data']

    anio = file[0:4]
    mes = file[4:6]
    dia = file[6:8]
    hora = file[8:12]
    fecha=anio+mes+dia+hora

    #el tercer parametro es la altura maxima hasta la cual 
    #voy a calcular las cosas (6000, aprox. 240km rango)
    #2000 es aprox. 120km rango)
    
    end_gate, start_ray, end_ray = ran.det_process_range(RADAR,0,2000) 
    dr = (RADAR.range['data'][1] - RADAR.range['data'][0]) / 1000.0

    at = np.zeros(maskphi.shape, dtype='float32')

    at=alpha*np.abs(kdp)
    z_cor=maskref+at

    #Enmascarar los datos invalidos con umbrales y con rho
    
    at[at<0]=np.nan
    at_ok = np.ma.masked_invalid(at)
    
  
    mask_at=np.ma.masked_where(np.isnan(RADARkdp.fields['mask_kdp']['data']),at_ok)
    aa=np.ma.filled(mask_at,fill_value=np.nan)
    bb = np.ma.masked_invalid(aa)
    
    mask_zcor=np.ma.masked_where(np.isnan(RADAR.fields['mask_ref']['data']),z_cor)
    dd=np.ma.filled(mask_zcor,fill_value=np.nan)
    ee = np.ma.masked_invalid(dd)

    #Agrego la nueva variable atenuacion con y sin mascara y la reflectividad
    #corregida por atenuacion

    RADAR.add_field_like('uncorrected_differential_phase','AT',at_ok)
    RADAR.add_field_like('uncorrected_differential_phase','mask_AT',bb)
    RADAR.add_field_like('uncorrected_differential_phase','ZCOR',z_cor)
    RADAR.add_field_like('uncorrected_differential_phase','mask_ZCOR',ee)

    
    ###### Grabo la variable atenuacion #####
    
    ATT=VariableUnica(RADAR, bb, 'dBZ','mask_att')
    
    ATT.add_field_like('mask_att','att',at_ok)
    
    ATT.add_field_like('mask_att','zcor',z_cor)
    
    ATT.add_field_like('mask_att','mask_zcor',ee)


    
    #***** Guardo el nuevo cfradial 
    pyart.io.write_cfradial(HOME + path + fecha+'att_simple.nc',ATT)

