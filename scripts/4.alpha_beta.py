# Este script calcula el ajuste lineal dado
# por una recta de regresion simple entre
# el zdr y Z con PHIDP
# para poder estimar el alpha y beta de la ecuacion 
# de atenuacion del metodo ZPHI

import numpy as np
import pyart
import matplotlib.pyplot as plt
import scipy.stats as stats
from os import listdir
from sys import argv
import glob

HOME='/home/romina/radar/doctorado/'

path=argv[1]


numero=len(glob.glob(HOME+path+"datos_PHIDP/*cor.nc"))

i=-1

alpha=np.zeros(numero)
beta=np.zeros(numero)
day=np.zeros(numero)

for fecha in list(set([f[0:12] for f in listdir(HOME + path+'/datos_PHIDP') if f.endswith('cor.nc') ])):
    
       
    i=i+1
    
    print(fecha)
    anio = fecha[0:4]
    mes = fecha[4:6]
    dia = fecha[6:8]
    hora = fecha[8:12]
    
    RADAR = pyart.io.read(HOME + path + 'datos_PHIDP/' + fecha + 'cor.nc')

    RADARkdp = pyart.io.read(HOME + path + 'datos_KDP/' + fecha + 'kdp.nc')

    RADARdelta = pyart.io.read(HOME + path + 'datos_DELTA/' + fecha + 'delta.nc')
    
    dbz= RADAR.fields['dBZ']['data']
    rho = RADAR.fields['RhoHV']['data']
    uphi = RADAR.fields['uncorrected_differential_phase']['data']
    phidp=RADAR.fields['maskPhi']['data']
    ref=RADAR.fields['mask_ref']['data']
    zdr = RADAR.fields['ZDR']['data']
    
    kdp=RADARkdp.fields['mask_kdp']['data']
    
    delta=RADARdelta.fields['delta']['data']

    
    ref[RADARdelta.fields['delta']['data']>5]=np.nan 
    ref[RADAR.fields['RhoHV']['data']<0.95]=np.nan
    ref[RADARkdp.fields['mask_kdp']['data']>2]=np.nan
    ref[RADARkdp.fields['mask_kdp']['data']<1]=np.nan
    
    phidp[RADARdelta.fields['delta']['data']>5]=np.nan 
    phidp[RADAR.fields['RhoHV']['data']<0.95]=np.nan
    phidp[RADARkdp.fields['mask_kdp']['data']>2]=np.nan
    phidp[RADARkdp.fields['mask_kdp']['data']<1]=np.nan
    
    zdr[RADARdelta.fields['delta']['data']>5]=np.nan 
    zdr[RADAR.fields['RhoHV']['data']<0.95]=np.nan
    zdr[RADARkdp.fields['mask_kdp']['data']>2]=np.nan
    zdr[RADARkdp.fields['mask_kdp']['data']<1]=np.nan
    zdr[RADAR.fields['ZDR']['data']<-50]=np.nan
    
    pru_ref=ref.flatten()
    pru_phi=phidp.flatten()
    pru_zdr=zdr.flatten()
    
    mask = ~np.isnan(pru_ref) & ~np.isnan(pru_phi)
    
    try:
        slope, intercept, r_value, p_value, std_err = stats.linregress(pru_phi[mask], pru_ref[mask])
    except:
        slope=np.nan
    
    alpha[i]=slope
    day[i]=fecha
    print(alpha)
           
    fig = plt.figure(figsize=(10,7))

    plt.plot(pru_phi[mask], pru_ref[mask], 'o', label='original data')
    plt.plot(pru_phi[mask], intercept + slope*pru_phi[mask], 'r', label='fitted line')
    plt.xlabel('Differential Propagation Phase [deg]')
    plt.ylabel('Horizontal Reflectivity [dBZ]')
    plt.title('Attenuation Correction Alpha '+dia+'-'+mes+'-'+anio+' '+hora+'UTC' ,fontsize=13)
    
    plt.savefig(HOME + path + fecha + "Z_PHI.png",dpi=150)
        
    plt.close(fig)
    
    mask_zdr = ~np.isnan(pru_zdr) & ~np.isnan(pru_phi)
    
    try:
        slope_b, intercept_b, r_value_b, p_value_b, std_err_b = stats.linregress(pru_phi[mask_zdr], pru_zdr[mask_zdr])
    except:
        slope_b=np.nan
        
    beta[i]=slope_b
    print(beta)
    
    fig=plt.figure(figsize=(10,7))
    plt.plot(pru_phi[mask_zdr], pru_zdr[mask_zdr], 'o', label='original data')
    plt.plot(pru_phi[mask_zdr], intercept_b + slope_b*pru_phi[mask_zdr], 'r', label='fitted line')
    plt.xlabel('Differential Propagation Phase [deg]')
    plt.ylabel('Differential Reflectivity [dB]')
    plt.title('Attenuation Correction Beta '+dia+'-'+mes+'-'+anio+' '+hora,fontsize=13)
    
    plt.savefig(HOME + path + fecha + "ZDR_PHI.png",dpi=150)
        
    plt.close(fig)

#print(day,alpha,beta)    
#np.savetxt(HOME+path+'alpha_prueba.out',(day,alpha,beta), delimiter=';')
