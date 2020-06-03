#Este script genera una mascara de 0 y 1 teniendo en cuenta el Rho
#y el ZDR. 
#Aplica la mascara al phidp y al dbz usando la primera y segunda elevacion
#y los desvios del phidp 

import matplotlib as mpl
import numpy as np
import pyart
import matplotlib.pyplot as plt
from os import listdir
from listado import * 
from sys import argv
import rango as ran
from funciones import promediomovil

##### LECTURA DE DATOS RAINBOW #########

HOME='/home/romina/radar/doctorado/'

path=argv[1]

for fecha in list(set([f[0:12] for f in listdir(HOME + path) if f.endswith('.vol') ])):
    uphidp, ref, cc, difref, vvel, ref_u = loadData(fecha,HOME,path)
    
    
    anio = uphidp.time['units'][14:18]
    mes = uphidp.time['units'][19:21]
    dia = uphidp.time['units'][22:24]
    hora = uphidp.time['units'][25:27] + uphidp.time['units'][28:30]
    fecha=anio+mes+dia+hora
    
    print(fecha)
    
    #Cuantos elementos hay ray y bins hay y cuales son las elevaciones
    nb1=ref.ngates
    nr1=ref.nrays
    levels=np.unique(ref.elevation['data'])
    
    dbz=ref.fields['reflectivity']['data']
    rho=cc.fields['cross_correlation_ratio']['data']
    zdr=difref.fields['differential_reflectivity']['data']

    
    mascara=np.zeros((nr1,nb1))
    dbz=np.ma.filled(dbz,fill_value=np.nan)
    rho=np.ma.filled(rho,fill_value=np.nan)
    zdr=np.ma.filled(zdr,fill_value=np.nan)

    
    try:
        for i in range(0,nr1):
            rho_h=rho[i,:]
    except:
        nr1=nr1-1
    
    for i in range(0,nr1):
        dbz_h=dbz[i,:]
        rho_h=rho[i,:]
        zdr_h=zdr[i,:]
        for j in range(0,nb1):
            if np.isnan(dbz_h[j]) or rho_h[j]<0.7 or np.isnan(rho_h[j]):
                mascara[i,j]=0
            else:
                if dbz_h[j]<20 and rho_h[j]<0.95 and zdr_h[j]>2:
                    mascara[i,j]=0
                else:
                    mascara[i,j]=1
    
    ### Agrego la variable mascara a la variable reflectividad            
    ref.add_field_like('reflectivity','mask',mascara,replace_existing='True')
    
    ### Tomo los datos de cada elevacion
                
    ref_ele=ref.extract_sweeps([0])
    ref_ele_ar=ref.extract_sweeps([1])
    uphi_ele=uphidp.extract_sweeps([0])
    uphi_ele_ar=uphidp.extract_sweeps([1])
    cc_ele=cc.extract_sweeps([0])
    difref_ele=difref.extract_sweeps([0])    
 
    mask_ab=ref_ele.fields['mask']['data']
    mask_ar=ref_ele_ar.fields['mask']['data']
    dbz_ab=ref_ele.fields['reflectivity']['data']
    uphi_ab=uphi_ele.fields['uncorrected_differential_phase']['data']
    rho_ab  = cc_ele.fields['cross_correlation_ratio']['data']
    zdr_ab = difref_ele.fields['differential_reflectivity']['data']

    #Cuantos elementos hay ray y bins hay
    nb=ref_ele.ngates
    nr=ref_ele.nrays 
    
    #Verifico en que posicion esta cada haz para las elevaciones 1 y 2

    ray_ar=np.zeros(nr)
    ray_ab=np.zeros(nr)

    for j in range(0,nr):
        desaz=j
        sweep=0

        i1_ab = ref_ele.sweep_start_ray_index['data'][sweep]
        i2_ab = ref_ele.sweep_end_ray_index['data'][sweep]
        azimuth_ab=ref_ele.azimuth['data'][i1_ab:i2_ab]
        azdiff_ab=abs(azimuth_ab-desaz)
        ray_num_ab = azdiff_ab.argmin()+i1_ab

        #En que posicion esta el azimuth que quiero (arriba)
        i1_ar = ref_ele_ar.sweep_start_ray_index['data'][sweep]
        i2_ar = ref_ele_ar.sweep_end_ray_index['data'][sweep]
        azimuth_ar=ref_ele_ar.azimuth['data'][i1_ar:i2_ar]
        azdiff_ar=abs(azimuth_ar-desaz)
        ray_num_ar = azdiff_ar.argmin()+i1_ar

        ray_ab[j]=ray_num_ab
        ray_ar[j]=ray_num_ar         
    
    ############# MASCARA PARA LA REFLECTIVIDAD ########################
        
    #Aplicar la mascara con los datos de la elevacion 0
    dbz_ab_filter=np.zeros((nr,nb))

    for i in range(0,nr):
        mask_ab_haz=mask_ab[i,:]
        dbz_ab_haz=dbz_ab[i,:]
        for j in range(0,nb):
            dbz_ab_filter[i,j]=mask_ab_haz[j]*dbz_ab_haz[j]
            if dbz_ab_filter[i,j] == 0:
                dbz_ab_filter[i,j]=np.nan
            else:
                dbz_ab_filter[i,j]=dbz_ab_haz[j]              
                    
    #Aplico la mascara a la elevacion 0 con datos de la elevacion 1
    dbz_filter=np.zeros((nr,nb))

    for i in range(0,nr):
        x=int(ray_ab[i])
        y=int(ray_ar[i])
        dbz_ab_haz=dbz_ab_filter[x,:]
        mask_ar_haz=mask_ar[y,:]
        for j in range(0,nb):
            dbz_filter[x,j]=mask_ar_haz[j]*dbz_ab_haz[j]
            if dbz_filter[x,j] == 0:
                dbz_filter[x,j]=np.nan
                
    #Corrijo los valores que dieron 0 porque se repite el azimuth en los datos
    for i in range(0,nr):
        for j in range(0,nb):
            if dbz_filter[i,j]==0:
                #print 'a'
                dbz_filter[i,:]=dbz_ab_filter[i,:]     
                
    ############# MASCARA PARA UPHIDP ########################

    #Aplicar la mascara con los datos de la elevacion 0
    uphi_ab_filter=np.zeros((nr,nb))

    for i in range(0,nr):
        mask_ab_haz=mask_ab[i,:]
        uphi_ab_haz=uphi_ab[i,:]
        for j in range(0,nb):
            uphi_ab_filter[i,j]=mask_ab_haz[j]*uphi_ab_haz[j]
            if uphi_ab_filter[i,j] == 0:
                uphi_ab_filter[i,j]=np.nan
            else:
                uphi_ab_filter[i,j]=uphi_ab_haz[j]
                
    #Aplico la mascara a la elevacion 0 con datos de la elevacion 1
    uphi_filter=np.zeros((nr,nb))

    for i in range(0,nr):
        x=int(ray_ab[i])
        y=int(ray_ar[i])
        uphi_ab_haz=uphi_ab_filter[x,:]
        mask_ar_haz=mask_ar[y,:]
        for j in range(0,nb):
            uphi_filter[x,j]=mask_ar_haz[j]*uphi_ab_haz[j]
            if uphi_filter[x,j] == 0:
                uphi_filter[x,j]=np.nan
                
    #Corrijo los valores que dieron 0 porque se repite el azimuth en los datos
    for i in range(0,nr):
        for j in range(0,nb):
            if uphi_filter[i,j]==0:
                uphi_filter[i,:]=uphi_ab_filter[i,:]                                
    
    end_gate, start_ray, end_ray = ran.det_process_range(ref,0,2000) 
    
    dbz_filter[:,end_gate+1:nb]=np.nan
    uphi_filter[:,end_gate+1:nb]=np.nan

    ########### MASCARA A LOS NAN #######
       
    aa_ref = np.ma.filled(dbz_filter,fill_value=np.nan)
    bb_ref = np.ma.masked_invalid(aa_ref)
    
    aa_uphi=np.ma.filled(uphi_filter,fill_value=np.nan)
    bb_uphi = np.ma.masked_invalid(aa_uphi)
    
    
###### ACA AGREGO LO DEL DESVIO DEL PHIDP ######
    
    phidp_haz=np.zeros(nb)
    des_uphi=np.zeros(nb)
    desvios=np.zeros((nr,nb))

    for j in range(0,nr):
        phidp_haz=bb_uphi[j,:]
        for k in range(0,nb):
            desvio_uphi=np.std(phidp_haz[k:k+2]) #Calculo desvios de uphi en 3 puntos
            desvios[j,k]=desvio_uphi
            
    phi_s=np.zeros((nr,nb)) 
    s1=np.zeros(nb)
    s2=np.zeros(nb)

    for j in range(0,nr):
        s1=desvios[j,:]
        s2=np.zeros(nb)
        for h in range (0,nb):
            s2=promediomovil(s1,3)
            phi_s[j,0:1]=np.nan
            phi_s[j,1:479]=s2            
    
    cc_desvios=np.ma.filled(phi_s,fill_value=np.nan)
    
    mask_des=np.zeros((nr,nb))

    for i in range(0,nr):
        bb_h=cc_desvios[i,:]
        for j in range(0,nb):
            if np.isnan(bb_h[j]):
                mask_des[i,j]=0
            else:
                mask_des[i,j]=1
                
    uphi_filter_des=np.zeros((nr,nb))
    ref_filter_des=np.zeros((nr,nb))

## Aplico esta mascara de desvios al phidp y a la reflectividad

    for i in range(0,nr):
        mask_haz=mask_des[i,:]
        uphi_haz=bb_uphi[i,:]
        for j in range(0,nb):
            uphi_filter_des[i,j]=mask_haz[j]*uphi_haz[j]
            if uphi_filter_des[i,j] == 0:
                uphi_filter_des[i,j]=np.nan
            else:
                uphi_filter_des[i,j]=uphi_haz[j]
    
    jj_uphi = np.ma.filled(uphi_filter_des,fill_value=np.nan)
    kk_uphi = np.ma.masked_invalid(jj_uphi)
    
    
    for i in range(0,nr):
        mask_haz=mask_des[i,:]
        ref_haz=bb_ref[i,:]
        for j in range(0,nb):
            ref_filter_des[i,j]=mask_haz[j]*ref_haz[j]
            if ref_filter_des[i,j] == 0:
                ref_filter_des[i,j]=np.nan
            else:
                ref_filter_des[i,j]=ref_haz[j]
    
    jj_ref = np.ma.filled(ref_filter_des,fill_value=np.nan)
    kk_ref = np.ma.masked_invalid(jj_ref)                            
    
    ######## Agrego las nuevas variables filtradas a una nueva y las viejas a un solo archivo nc #####
    
    uphi_ele.add_field_like('uncorrected_differential_phase','mask_uphi',kk_uphi,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','dBZ',dbz_ab,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','RhoHV',rho_ab,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','ZDR',zdr_ab,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','mask_ref',kk_ref,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','mascara_ab',mask_ab,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','mascara_ar',mask_ar,replace_existing='True')
    uphi_ele.add_field_like('uncorrected_differential_phase','mascara_des',mask_des,replace_existing='True')
    

    
    uphi_ele.add_field_like('uncorrected_differential_phase','desvios',cc_desvios,replace_existing='True')

    
    ####### GUARDO EL NUEVO NETCDF ############ 
     
    pyart.io.write_cfradial(HOME + path + fecha+'.nc',uphi_ele) 
 
 
    
########################### GRAFICOS DE FILTER UPHI Y FILTER DBZ ####################################
    
    ele=str(levels[sweep])

    display_uphi = pyart.graph.RadarDisplay(uphi_ele)

    fig = plt.figure(figsize=(9,7))


    xrange  = [-120,120]
    yrange  = [-120,120]
    anillos = [30,60,90,120]
    Rmax=120.0
        
    xlabel = 'Distance (km)'
    ylabel = 'Distance (km)'

    display_uphi.plot_ppi('mask_ref',colorbar_label='dBZ',cmap='pyart_NWSRef',axislabels=(xlabel,ylabel),vmin=0,vmax=70)
    plt.title('Horizontal Reflectivity (filtered) ' +dia+'-'+mes+'-'+anio+' '+hora+' '+ele+' deg',fontsize=15,y=1.02)
    display_uphi.plot_range_rings(anillos,lw='0.7')
    display_uphi.set_limits(xrange,yrange)
    display_uphi.plot_cross_hair(1)

    plt.tight_layout()

    plt.savefig(HOME + path + fecha + "filter_dbz.png",dpi=150)
        
    plt.close(fig)
    
    fig = plt.figure(figsize=(9,7))

    xrange  = [-120,120]
    yrange  = [-120,120]
    anillos = [30,60,90,120]
    Rmax=120.0
        
    xlabel = 'Distance (km)'
    ylabel = 'Distance (km)'

    display_uphi.plot_ppi('mask_uphi',colorbar_label='uPhiDP',cmap='pyart_NWSRef',axislabels=(xlabel,ylabel),vmin=0,vmax=360)
    plt.title('Uncorrected Differential Phase (filtered) ' +dia+'-'+mes+'-'+anio+' '+hora+' '+ele+' deg',fontsize=15,y=1.02)
    display_uphi.plot_range_rings(anillos,lw='0.7')
    display_uphi.set_limits(xrange,yrange)
    display_uphi.plot_cross_hair(1)

    plt.tight_layout()
    
    plt.savefig(HOME + path + fecha + "filter_uphi.png",dpi=150)
        
    plt.close(fig)
                    
