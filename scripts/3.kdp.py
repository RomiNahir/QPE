# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:52:45 2017

@author: romina
"""

import numpy as np
import pyart
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob, os
from sys import argv
from funciones import VariableUnica


HOME='/media/romina/datos/ScriptsTesis/Doctorado'

path=argv[1]


os.chdir(HOME+path)
for file in glob.glob("*cor.nc"):
    print(file)
    
    RADAR = pyart.io.read(HOME + path + file)
    
    dbz= RADAR.fields['dBZ']['data']
    maskphi=RADAR.fields['maskPhi']['data']
    
    anio = file[0:4]
    mes = file[4:6]
    dia = file[6:8]
    hora = file[8:12]
    fecha=anio+mes+dia+hora

    nb=RADAR.ngates
    nr=RADAR.nrays
    
    kdp=np.zeros((nr,nb))

    for j in range(0,nr):
        phi_kdp=maskphi[j,:]
    
        for i in range(0,nb): 
            s1=max(0,i-2)
            s2=min(i+2,nb-1)
            r=[x*1. for x in range(s2-s1+1)]
            x=2.*0.5*np.array(r)
            y=phi_kdp[s1:s2+1]
            a=len(y[np.where(y>0)])
            b=np.std(y)
  
            if a==5:# and b<20:
                ajuste=np.polyfit(x,y,1)
                kdp[j,i]=ajuste[0]
            else:
                kdp[j,i]=np.nan

    #Enmascarar los datos invalidos con umbrales y con la mascara

#    kdp[kdp>15]=np.nan
#    kdp[kdp<-10]=np.nan

    kdp_ok = np.ma.masked_invalid(kdp)
       
#    mask_kdp=np.ma.masked_where(np.isnan(RADAR.fields['mask_uphi']['data']),kdp_ok)
    mask_kdp=np.ma.masked_where(np.isnan(RADAR.fields['maskPhi']['data']),kdp_ok)
	
    aa=np.ma.filled(mask_kdp,fill_value=np.nan)

    bb = np.ma.masked_invalid(aa)
        
    #Agrego la nueva variable kdp

 
    KDP=VariableUnica(RADAR, bb, 'dBZ','mask_kdp')
    
    KDP.add_field_like('mask_kdp','kdp',bb)

    
    #***** Guardo el nuevo cfradial 
    pyart.io.write_cfradial(HOME + path + fecha+'kdp.nc',KDP)                   
    
    
    #######################  GRAFICO   ##################
    
    font = {'family' : 'sans-serif',
        'size'   : 19}

    mpl.rc('font', **font)
    
    display_kdp = pyart.graph.RadarDisplay(RADAR)

    xrange  = [-120,120]
    yrange  = [-120,120]
    anillos = [30,60,90,120]
    Rmax=120.0

    xlabel = 'Distance (km)'
    ylabel = 'Distance (km)'
       
    fig = plt.figure(figsize=(9,7))

    display_kdp.plot_ppi('mask_kdp',colorbar_label='deg/km',cmap='pyart_Theodore16',axislabels=(xlabel,ylabel),vmin=-5,vmax=10)
    plt.title('$K_{DP}$ ' +dia+'-'+mes+'-'+anio+' '+hora+' '+'0.5 deg',fontsize=19,y=1.02)
    display_kdp.plot_range_rings(anillos,lw='0.7')
    display_kdp.set_limits(xrange,yrange)
    display_kdp.plot_cross_hair(1)

    plt.tight_layout()
    
    plt.savefig(HOME + path + fecha + "kdp.png",dpi=150)
        
    plt.close(fig)
