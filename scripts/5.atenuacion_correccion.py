# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:52:45 2017

@author: romina
"""


### Este script calcula la atenuacion y la correcion en reflectividad
### Es importante que previamente sepa cual es el alpha
### que se va a usar porque puede cambiar según la situación
### el alpha se calcula en el paso 4.alpha_beta.py


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

HOME='/media/romina/datos/Doctorado/'

# el path es donde estan los corregidos de phidp

path=argv[1]


os.chdir(HOME+path)
for file in glob.glob("*cor.nc"):
    print(file)
    
    RADAR = pyart.io.read(HOME + path + file)
    
    levels=np.unique(RADAR.elevation['data'])

    #el coeficiente alpha para banda C calculados con el alpha_beta
    #coeficiente b, que varia entre 0.6 y 0.9 segun ryzkhov 2014 y es 0.99 segun bringi libro para banda C 

    b=0.8    #coeficiente b
    alfa=0.0654  #coeficiente alfa
    
    rho = RADAR.fields['RhoHV']['data']
    maskphi=RADAR.fields['maskPhi']['data']
    maskref=RADAR.fields['mask_ref']['data']


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
    atten = np.zeros(maskphi.shape, dtype='float32')

    rhv_nan=0.6
    phidp_nan=0.0

    # Aplico una mascara en los datos dependiendo de rho
    is_cor = rho > rhv_nan
    is_phi = maskphi > phidp_nan
    is_good = np.logical_and(is_cor, is_phi)

    for i in range(start_ray, end_ray-1):
        
        ray_phi = maskphi[i, 0:end_gate]
        ray_dbz = maskref[i, 0:end_gate]
        
        ultimos = np.where(is_good[i, 0:end_gate])[0][-5:]
        primeros = np.where(is_good[i,0:50])[0][-10:]
        
        last_six_good = np.where(is_good[i, 0:end_gate])[0][-6:]

        try:
            dif=np.nanmedian(ray_phi[ultimos])-np.nanmin(ray_phi[primeros])
            
        except:
            dif = np.median(ray_phi[last_six_good])
        
        if dif>0:
            PIA = alfa * dif
        else:
            PIA= alfa * np.median(ray_phi[last_six_good])

        #promedio movil con convolucion
        sm_refl = smooth_and_trim(ray_dbz, window_len=5) 

        #tranformo dbz a mm6 m-3 y ya lo multiplico por la constante b
        reflectividad_lineal = 10.0 ** (0.1 * b * sm_refl) 
        
        a=(np.exp(0.23*b*PIA))-1

        I_all = cumtrapz(0.46 * b * dr * reflectividad_lineal[:: -1]) 
            
        I_all = np.append(I_all, I_all[-1])[::-1]
        
        # Calculo de la atenuacion especifica
        at[i, 0:end_gate] = (reflectividad_lineal * a /(I_all[0] + a * I_all))
        
          
        # Correccion de Z por atenuacion
        atten[i, :-1] = cumtrapz(at[i, :]) * dr * 2.0
        atten[i, -1] = atten[i, -2]
   
        z_cor=maskref+atten

    #Enmascarar los datos invalidos con umbrales y con rho
    
    at[at<0]=np.nan
    at_ok = np.ma.masked_invalid(at)
    
  
    mask_at=np.ma.masked_where(np.isnan(RADAR.fields['mask_uphi']['data']),at_ok)
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
    pyart.io.write_cfradial(HOME + path + fecha+'att.nc',ATT)

    ################### GRAFICOS ###################
    sweep=0
    ele=str(levels[sweep])

    font = {'family' : 'sans-serif',
        'size'   : 19}

    mpl.rc('font', **font)

    display = pyart.graph.RadarDisplay(RADAR)

    xrange  = [-120,120]
    yrange  = [-120,120]
    anillos = [30,60,90,120]
    Rmax=120.0

    xlabel = 'Distance (km)'
    ylabel = 'Distance (km)'

    #### SOLO ATENUACION #####
    fig = plt.figure(figsize=(9,7))
    
    display.plot_ppi('mask_AT',colorbar_label='dB/km',cmap='pyart_Carbone17',axislabels=(xlabel,ylabel),vmin=0,vmax=1)
    plt.title('Specific Attenuation ($A_{H}$) ' +dia+'-'+mes+'-'+anio+' '+hora+' '+ele+' deg',fontsize=19,y=1.02)
    display.plot_range_rings(anillos,lw='0.7')
    display.set_limits(xrange,yrange)
    display.plot_cross_hair(1)

    plt.tight_layout()
    
    plt.savefig(HOME + path + fecha + ele+"att.png",dpi=150)
        
    plt.close(fig)    
    
    #### SOLO REFLECTIVIDAD CORREGIDA ####
    fig = plt.figure(figsize=(9,7))

    display.plot_ppi('ZCOR',colorbar_label='dBZ',cmap='pyart_NWSRef',axislabels=(xlabel,ylabel),vmin=0,vmax=70)
    plt.title('$Z_{H}$ (corrected) ' +dia+'-'+mes+'-'+anio+' '+hora+' '+ele+' deg',fontsize=19,y=1.02)
    display.plot_range_rings(anillos,lw='0.7')
    display.set_limits(xrange,yrange)
    display.plot_cross_hair(1)

    plt.tight_layout()
    
    plt.savefig(HOME + path + fecha + ele+"dbz_cor.png",dpi=150)
        
    plt.close(fig)
