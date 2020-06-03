import glob
import re
import pyart



def getfilename(pattern):
    try:
        if pattern.endswith('V.vol') and not pattern.endswith('RhoHV.vol'):
            return [p for p in glob.glob(pattern) if pattern.endswith('V.vol') and not pattern.endswith('RhoHV.vol')][0]
        else:
            return glob.glob(pattern)[0]
    except Exception as e:
        #print "Error:",e
        return None

def loadData(fecha_hora,HOME,path):
    full_path = HOME + path
    try:
        uphidp = pyart.aux_io.read_rainbow_wrl(getfilename( "%s%s%s*uPhiDP.vol" % (HOME, path, fecha_hora)))
    except:
        uphidp = None
        
    try:    
        ref = pyart.aux_io.read_rainbow_wrl(getfilename( "%s%s%s*dBZ.vol" % (HOME, path, fecha_hora)))
    except:
        ref = None
        
    try:    
        cc = pyart.aux_io.read_rainbow_wrl(getfilename( "%s%s%s*RhoHV.vol" % (HOME, path, fecha_hora)))
    except:
        cc = None
        
    try:    
        difref = pyart.aux_io.read_rainbow_wrl(getfilename( "%s%s%s*ZDR.vol" % (HOME, path, fecha_hora)))
    except:
        difref = None
        
    try:    
        vvel = pyart.aux_io.read_rainbow_wrl(getfilename( "%s%s%s*[0-9]V.vol" % (HOME, path, fecha_hora)))
    except:
        vvel = None
        
    try:    
        ref_u = pyart.aux_io.read_rainbow_wrl(getfilename( "%s%s%s*dBuZ.vol" % (HOME, path, fecha_hora)))
    except:
        ref_u = None        
    
    return uphidp, ref, cc, difref, vvel,ref_u 
