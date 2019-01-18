
"""
GENERATES POLARIZED SYNCHROTRON FOREGROUNDS AT DIFFERENT FREQUENCIES

"""

import numpy as np
import healpy as hp
import os,sys
import random
from astropy.io import fits
import scipy.optimize 
from scipy import stats
import invRM as iRM 
import genRMmap as gRM

c = 299792458.0 #speed of light



def arcmin2rad(arcmin):
    return arcmin/60/180*np.pi


def jyb2Tb(jyb, res,freq, rmtf):
    "convert Jy beam^-1 rmtf-1 in Tb[K] given res in arcmin and freq in MHz"
    rad=arcmin2rad(res)
    #omega=np.pi*(rad)**2/2/np.log(2)
    omega=rad**2
    return (jyb/omega*rmtf)*1e-26*(3*1e8)**2/(2*1.38*1e-23*(freq*1e6)**2)


  
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] DATA_MWA')
    o.set_description(__doc__)
    o.add_option('--seed',dest='seed',default=1001,help='seed for generating RM maps')
    o.add_option('--start',dest='start',default=50,help='Starting frequency in MHz')
    o.add_option('--stop',dest='stop',default=200,help='Stopping frequency in MHz')
    o.add_option('--nchan',dest='nchan',default=150,help='Number of frequency channels')
    o.add_option('--nside',dest='nside',default=128,help='Resolution of the map')
    o.add_option('--filepath',dest='filepath',help='Path of Output map fits files')
    o.add_option('--storeRM',dest='storeRM',action='store_true', help='Set option to store RM space Q,U maps')
    o.add_option('--readRM',dest='readRM',action='store_true', help='Set option to use RM maps from fits file')
    o.add_option('--fileq',dest='fileq',default='mapRM_Q',help='file RM Q maps from fits DEFAULT:mapRM_Q_nside[nside]_[RM].fits')
    o.add_option('--fileu',dest='fileu',default='mapRM_U',help='file RM U maps from fits DEFAULT:mapRM_U_nside[nside]_[RM].fits')

    opts, args = o.parse_args(sys.argv[1:]) 

    print opts, args 

    nside = int(opts.nside)
    npix=hp.nside2npix(nside)
    seed=int(opts.seed)
    filepath=str(opts.filepath)
    start,stop,nchan = np.float(opts.start),np.float(opts.stop),np.float(opts.nchan)
    RM= np.load(args[0])['RM']     
    sigma = np.load(args[0])['sigma']  
    alpha = np.load(args[0])['alpha']  
    frequency = np.linspace(start,stop,nchan)
    np.savez(open(filepath+'input_info.npz','wb'),frequency=frequency,nside=nside,seed=seed)
    
    ##to convert in K
    convert=jyb2Tb(1,15.6,189, 4.3)  ##See Bernardi et al 2013
    
    print "Input read"
    

    if opts.readRM:
        print "Reading RM maps from input files"
        RM_cube=np.zeros((len(RM),npix,2),dtype=float)
        for rm in range(len(RM)):
            fileq=str(opts.fileq)+'_nside'+str(nside)+'_seed'+str(seed)+'_'+str(rm)+'.fits' #Note input file need to have nside=NSIDE
            fileu=str(opts.fileu)+'_nside'+str(nside)+'_seed'+str(seed)+'_'+str(rm)+'.fits' 
            RM_cube[rm,:,0]=hp.read_map(fileq)
            RM_cube[rm,:,1]=hp.read_map(fileu)

            print "Q,U maps read for RM=", RM[rm] 



    else:
        print "Generating RM full sky maps.."
        ##setting the seed
        ## Here just used a simple generation for seed_vec 
        seed_vec=np.arange(seed,seed+len(RM),1)
        RM_cube=np.zeros((len(RM),npix,2),dtype=float)

        for rm in range(len(RM)):
            map_QU=np.zeros((nside,2))
            map_QU=gRM.genRMmapQU(alpha[rm], nside,seed_vec[rm],sigma[rm]) 
            RM_cube[rm,:,:]=map_QU
            print "Q,U maps genereted for RM=", RM[rm] 

    	    if opts.storeRM: 

                fname_Q=filepath+'mapRM_Q_nside'+str(nside)+'_seed'+str(seed)+'_'+str(rm)+'.fits'
                hp.fitsfunc.write_map(fname_Q, map_QU[:,0])
                fname_U=filepath+'mapRM_U_nside'+str(nside)+'_seed'+str(seed)+'_'+str(rm)+'.fits'
                hp.fitsfunc.write_map(fname_U, map_QU[:,1])      
                ## Note all file are stored with same seed value (the input one)

   
            print "RM full sky maps generated."
   
    
    print "Apply inverse RM"
    
 
    
    for f in frequency:
        p_map=iRM.P_l2(npix, RM_cube[:,:,0], RM_cube[:,:,1], RM, f)
        print "Computed freq=",f

   
        map_Q=p_map.real*convert
        map_U=p_map.imag*convert

        freq,df=str(f).split('.')

        fname_Q=filepath+'map_Q_nside'+str(nside)+'_seed'+str(seed)+'_freq'+freq+'p'+df[:3]+'.fits'
        hp.fitsfunc.write_map(fname_Q, map_Q)
        fname_U=filepath+'map_U_nside'+str(nside)+'_seed'+str(seed)+'_freq'+freq+'p'+df[:3]+'.fits'
        hp.fitsfunc.write_map(fname_U, map_U)
      
