# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:52:45 2017

@author: romina
"""
# Lee la salida de los datos netcdf emascarados con "1.mascara"
# Corrige el phidp, guarda todas las variable de la elevacion mas
# mas baja en un netcdf y luego hace un grafico hasta 120km
# toma umbrales de desvios y percentiles, hace unfolding
# suaviza la variable
# calcula el phi del sistema y lo resta
# corrige por radiales y le aplica una mascara 
# tomando umbrales de reflectividad y coeficiente de correlacion

import matplotlib as mpl
#mpl.use('agg')
import numpy as np
import pyart
import matplotlib.pyplot as plt
from funciones import promediomovil
import sistema as sis
import glob, os
from sys import argv
import rango as ran

HOME='/home/romina/radar/doctorado'


path=argv[1]

os.chdir(HOME+path)
for file in glob.glob("*.nc"):
    print(file)
    
    RADAR = pyart.io.read(HOME + path + file)
    
    zdr = RADAR.fields['ZDR']['data']
    dbz= RADAR.fields['mask_ref']['data']
    rho = RADAR.fields['RhoHV']['data']
    uphi = RADAR.fields['mask_uphi']['data']
    raw_dbz = RADAR.fields['dBZ']['data']
    
    nb=RADAR.ngates
    nr=RADAR.nrays
    
    end_gate, start_ray, end_ray = ran.det_process_range(RADAR,0,2000) 
    
    raw_dbz[:,end_gate+1:nb]=np.nan
    rho[:,end_gate+1:nb]=np.nan
    
    raw_dbz=np.ma.filled(raw_dbz,fill_value=np.nan)
    rho=np.ma.filled(rho,fill_value=np.nan)
    
    raw_dbz=np.ma.masked_invalid(raw_dbz)
    rho=np.ma.masked_invalid(rho)
    
    anio = file[0:4]
    mes = file[4:6]
    dia = file[6:8]
    hora = file[8:12]
    fecha=anio+mes+dia+hora
    
    
    
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
    
    ######## Suaviza el phi unfolding y chequeado
    ######## Suavizado con una media movil de x puntos (metodo convolucion)
    
    phi_s=np.zeros((nr,nb)) 
    s1=np.zeros(nb)
    s2=np.zeros(nb)
    
    for j in range(0,nr):
        s1=phi_cor[j,:]
        for h in range (0,nb):
            s2=promediomovil(s1,5)
        phi_s[j,0:1]=np.nan
        phi_s[j,479]=np.nan
        phi_s[j,2:478]=s2  
    
    ############# Calculo el phi del sistema y se lo resta
    
    system_zero = sis.det_sys_phase(RADAR,ncp_lev=0.6,rhohv_lev=0.6,
                                    ncp_field='RhoHV', rhv_field='RhoHV'
                                    , phidp_field='uncorrected_differential_phase')
    
    #**************************************************************        
    # Aca resto el phi del sistema al phi corregido y suavizado
    
    phi_final=np.zeros((nr,nb))
    phi_err=np.zeros((nr,nb))
    
    print(system_zero)
    
    
    for k in range(0,nr):
        for p in range(0,nb):
            phi_err[k,p]=np.nan
    
    try:
        for j in range(0,nr):
            for i in range(0,nb):
                phi_final[j,i]=phi_s[j,i]-system_zero
    except:
        phi_final=phi_err
    
    #### Agrego la nueva variable phidp corregida por folding y phi del sistema
    
    RADAR.add_field_like('mask_uphi','PhiDPcor',phi_final,replace_existing='True')
    
    #Correccion por radiales imagen completa
    
    des_rad=np.zeros((nr,nb))
    phi_rad=np.zeros((nr,nb))
    d2=np.zeros(nr)
    d3=np.zeros(nr)
    
    for k in range(0,nb):
        d2=phi_final[:,k]
        for j in range (0,nr-1):
            desvio_radial=np.std(d2[j-1:j+1]) 
            dif1=abs(d2[j]-d2[j-1])
            dif2=abs(d2[j]-d2[j+1])
            if desvio_radial > 20 and dif1 > 100 and dif2 > 100:
                f=[d2[j-1],d2[j+1]]
                d3[j]=np.mean(f)
            else:
                d3[j]=d2[j]
        phi_rad[:,k]=d3
        
    #Corrijo los valores que dieron 0 porque se repite el azimuth en los datos
    for i in range(0,nr):
        for j in range(0,nb):
            if phi_rad[i,j]==0:
                phi_rad[i,:]=phi_final[i,:]      
        
    
    RADAR.add_field_like('mask_uphi','PhiDPrad',phi_rad,replace_existing='True')
    
    #enmascaro los nan del phidprad y agrego el phidp final corregido en forma radial
    
    mask_phi_rad = np.ma.masked_invalid(phi_rad)
    
    RADAR.add_field_like('mask_uphi','maskPhiDPrad',mask_phi_rad,replace_existing='True') 
    
    std=np.zeros(nr)
    final=np.zeros((nr,nb))
    
    for i in range(0,nr):
        haz=mask_phi_rad[i,2:nb]
        try:
            std[i]=np.nanstd(haz)
        except:
            std[i]=0
    
    mask_phi=np.ma.masked_where(np.isnan(RADAR.fields['mask_uphi']['data']),mask_phi_rad)
    aa=np.ma.filled(mask_phi,fill_value=np.nan)
    bb = np.ma.masked_invalid(aa)
    
    #para que no se borren los pixels por detras con datos uniformes
    #se analiza esto, el desvio standard de cada rayo
    #asi no se borran con la mascara
    
    for i in range(0,nr):   
        if std[i]<30:
            final[i,:]=bb[i,:]
        else:
            final[i,:]=mask_phi_rad[i,:]
                            
    final=np.ma.masked_invalid(final)
    
    final[:,end_gate+1:nb]=np.nan
    final=np.ma.masked_invalid(final)
    
    minimo=np.abs(np.nanmin(final[:,2:]))
    print(minimo)
    
    if minimo!=system_zero:
        phi=final+minimo
    else:
        phi=final
    #phi[:,end_gate+1:nb]=np.nan
    #phi=np.ma.masked_invalid(phi)
    
    RADAR.add_field_like('mask_uphi','maskPhi',phi,replace_existing='True')
    
    #***** Guardo el nuevo cfradial 
    
    pyart.io.write_cfradial(HOME + path + fecha+'cor.nc',RADAR)
    
