# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:52:45 2017

@author: romina
"""
# Calcula el diferencial de cambio de fase
# de retrodispersion (backscatter diferential phase) con el phidp 
# enmascarado, utilizando  
# una media de 8 puntos 
# moviles

# la carpeta de path es donde estan los ENMASCARADO


import numpy as np
import rango as ran
import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
from funciones import promediomovil
from funciones import VariableUnica
import glob, os
from sys import argv

HOME='/home/romina/radar/doctorado/'

path=argv[1]

os.chdir(HOME+path)
for file in glob.glob("*.nc"):
    print(file)
    
    RADAR = pyart.io.read(HOME + path + file)

    uphi = RADAR.fields['mask_uphi']['data']

    anio = file[0:4]
    mes = file[4:6]
    dia = file[6:8]
    hora = file[8:12]
    fecha=anio+mes+dia+hora

    nb=RADAR.ngates
    nr=RADAR.nrays
    
    end_gate, start_ray, end_ray = ran.det_process_range(RADAR,0,2000) 
    
    #### Calculo los desvios y las media para cada haz  ######
    
    phidp_haz=np.zeros(nb)
    des_uphi=np.zeros(nb)
    desvios=np.zeros((nr,nb))
    medias=np.zeros((nr,nb))
    
    for j in range(0,nr):
        phidp_haz=uphi[j,:]
    
        for k in range(0,nb):
            desvio_uphi=np.std(phidp_haz[k:k+2]) #Calculo desvios de uphi en 3 puntos
    
            des_uphi[k]=desvio_uphi
    
            media_desvio=np.nanmean(des_uphi[k-20:k])  #calculo media de esos desvios en 20 puntos 
    
            desvios[j,k]=desvio_uphi
            medias[j,k]=media_desvio
        
    ##### Calculo los percentiles para cada haz #########
    
    per_m=np.zeros(nr)
    per_d=np.zeros(nr)
    
    for h in range(0,nr):
        desvio_haz=desvios[h,:]
        media_haz=medias[h,:]
    
        per_d[h]=np.nanpercentile(desvio_haz,99) # np.arange(0, 100, 10))
        per_m[h]=np.nanpercentile(media_haz,90) #np.arange(0, 100, 10))
    
    #Comparo los desvios con los percentiles para poder determinar que datos tienen
    #errores y donde no se hará el unfolding
    
    phi_d=np.zeros((nr,nb))
    
    phidp_haz=np.zeros(nb)
    
    for j in range(0,nr):
        phidp_haz=uphi[j,:]
    
        for k in range(0,nb):
    
            if desvios[j,k] > per_d[j] and medias[j,k] > per_m[j]:
                phi_d[j,k]=np.nan
            else:
                phi_d[j,k]=phidp_haz[k]        
    
    ##### Aca hace el unfolding  ###########
    
    phi_cor=np.zeros((nr,nb)) #Asigno cero a la nueva variable phidp corregida
    
    v1=np.zeros(nb)  #Vector v1
    #v2=np.zeros(nb)  #Vector v2
    
    diferencia=280 # Valor que toma la diferencia entre uno y otro pixel dentro de un mismo azimuth
    
    for m in range(0,nr):
        v1=phi_d[m,:]
        v2=np.zeros(nb)
        for l in range(0,nb):
            a=v2[l-1]-v1[l]
            if np.isnan(a):
                v2[l]=v2[l-1]
            else:    
                if a>diferencia:
                    v2[l]=v1[l]+360
                    if v2[l-1]-v2[l]>100:   #esto es por el doble folding cuando es mayor a 700
                        v2[l]=v1[l]+v2[l-1]
                else:
                    v2[l]=v1[l]
        phi_cor[m,:]=v2
    
    phi_cor[phi_cor==0.0]=np.nan
    
    
    ######## Suavizado con una media movil de x puntos (metodo convolucion)

    phi_s=np.zeros((nr,nb)) 
    s1=np.zeros(nb)
    s2=np.zeros(nb)

    for j in range(0,nr):
        s1=phi_cor[j,:]
        s2=np.zeros(nb)
        for h in range (0,nb):
            s2=promediomovil(s1,8)
        phi_s[j,0:4]=np.nan
        phi_s[j,479]=np.nan
        phi_s[j,5:478]=s2
            
    #############
    # Backscatter Differential Phase.... Difference between raw PHIDP and
    # Smooth PHIDP (delta)  

    delta_1=phi_cor-phi_s
    
    ### Suavizado del delta  En este suavizado 
    ### se pierden algunos valores.
    ### Tener en cuenta esto si se utiliza.
    
    delta_s=np.zeros((nr,nb)) 
    d1=np.zeros(nb)
    d2=np.zeros(nb)

    for j in range(0,nr):
        d1=delta_1[j,:]
        d2=np.zeros(nb)
        for h in range (0,nb):
            d2=promediomovil(d1,5)
        delta_s[j,0:1]=np.nan
        delta_s[j,479]=np.nan
        delta_s[j,2:478]=d2
    
    
    #Enmascarar los datos invalidos con umbrales y con rho
    #delta=np.abs(delta)
    delta_s[np.abs(delta_s)>100]=np.nan
    delta_s[np.abs(delta_s)==0]=np.nan

    #percentil=np.nanpercentile(delta,90)
    #print(percentil)
    
    #delta[delta>percentil]=np.nan
    
    delta_s[:,end_gate+1:nb]=np.nan
    
    delta_ok = np.ma.masked_invalid(delta_s)
    print(np.min(delta_ok))


##### Agrego la variable backscatter ########

    RADAR.add_field_like('mask_uphi','DELTA',delta_ok)
    
    #Agrego la nueva variable kdp

 
    DELTA=VariableUnica(RADAR, delta_ok, 'dBZ','delta')
    
    DELTA.add_field_like('delta','backscatter',delta_ok)

    
    #***** Guardo el nuevo cfradial 
    pyart.io.write_cfradial(HOME + path + fecha+'delta.nc',DELTA)
