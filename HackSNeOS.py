# coding: utf-8
from __future__ import absolute_import

import os
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt

import sncosmo

import gedankenLSST
from LSSTmetrics import PerSNMetric
from lsst.sims.photUtils import BandpassDict
from analyzeSN import analyzelcFits as anf
from astropy.units import Unit


lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()

# sncosmo Bandpasses required for fitting
throughputsdir = os.getenv('THROUGHPUTS_DIR')

bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')

for band in bandPassList:

    # setup sncosmo bandpasses
    bandfname = banddir + "/total_" + band + '.dat'


    # register the LSST bands to the SNCosmo registry
    # Not needed for LSST, but useful to compare independent codes
    # Usually the next two lines can be merged,
    # but there is an astropy bug currently which affects only OSX.
    numpyband = np.loadtxt(bandfname)
    #print band
    sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                   trans=numpyband[:, 1],
                                   wave_unit=Unit('nm'),
                                   name=band)
    sncosmo.registry.register(sncosmoband, force=True)


# In[3]:

lsstCadence = deepcopy(gedankenLSST.LSSTReq)


# In[4]:

# function that produces expected value of variances for many shifted lightcurves
# sncosmo_lc


def shift_loop(delta_t0, lsst_obs):
    var = []

    for i in delta_t0:
        try:
            sn = PerSNMetric(summarydf=lsst_obs.summary,t0=i, raCol='ra', decCol='dec', lsst_bp=lsst_bp)
            var.append(1.0/sn.qualityMetric(Disp=1.0))
            #print 'variance!', var[i]
        except:
            print('I failed!')
    
    return np.mean(np.array(var))

def shift_loop_mcmc(delta_t0, lsst_obs):
    model = sncosmo.Model(source='salt2-extended')
    var = []

    for i in delta_t0:
        try:
            sn = PerSNMetric(summarydf=lsst_obs.summary,t0=i, raCol='ra', decCol='dec', lsst_bp=lsst_bp)
            data = sn.SNCosmoLC() 
            mcmc_out = sncosmo.mcmc_lc(data,model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(0.3, 0.7)})
            t = anf.ResChar.fromSNCosmoRes(mcmc_out)
            print(i,t.salt_samples().mu.std())
            var.append(t.salt_samples().mu.std()*t.salt_samples().mu.std())
        except:
            print('I failed!')
    
    return np.mean(np.array(var))

def cadence_loop(bumps,delta_t0):

    #delta_t0 = np.linspace(snLSST.SN.mintime(), snLSST.SN.maxtime(), 10)
    #print delta_t0
    #print bumps

    length = []
    mu_variance_per_bump = []
    for i in bumps:
        print i
        lsstCadence['bF'] = i
        
        lsst_obs = gedankenLSST.GSN_Obs(mjd_center=49570., 
                                    lsstrequirements=lsstCadence,
                                    ra=58., dec=-27.,
                                    timeWindow=[-130., 150.]) 
        length.append(len(lsst_obs.summary))
        mu_variance_per_bump.append(shift_loop(delta_t0, lsst_obs))
    
    return mu_variance_per_bump, length

def cadence_loop_mcmc(bumps,delta_t0):
    ## Could move this function into the previous function and have mcmc be an option.
    length = []
    mu_variance_per_bump = []
    for i in bumps:
        print 'In i loop, ', i
        exit()
        lsstCadence['bF'] = i
        
        lsst_obs = gedankenLSST.GSN_Obs(mjd_center=49570., 
                                    lsstrequirements=lsstCadence,
                                    ra=58., dec=-27.,
                                    timeWindow=[-130., 150.]) 
        length.append(len(lsst_obs.summary))
        mu_variance_per_bump.append(shift_loop_mcmc(delta_t0, lsst_obs)) ### This line is the only difference between the two.
    
    return mu_variance_per_bump, length


if __name__ == '__main__':


    ### CONTROL HERE: TIME STEP
    t0 = 49570.
    dt = 50.
    delta_t0 = np.linspace(t0-dt,t0+dt,10)
    print(delta_t0)
    #49540.0
    #49645.0
    
    #normalize Max-Min/ time window
    
    ## Using Maximum likelihood  methods
    #bumps = np.arange(0.1, 3, 0.2)
    #mvpb, length = cadence_loop(bumps,delta_t0)
    
    
    
    
    ### CONTROL HERE: NUMBER AND RANGE OF BUMPS
    bumps1 = np.arange(1.0, 2.0, 0.5)
    mvpb1, length1 = cadence_loop_mcmc(bumps1,delta_t0)
    
    np.savetxt('mu_variance_length.txt', zip(mvpb1,length1))
    
    
    
    normalize =(49645.0 - 49540.0)/(130.+ 150.)
    
    x = np.linspace(0.1, 80, 1000)
    
    plt.figure(figsize = (8,6))
    plt.plot(np.array(length1)*normalize, mvpb1, 'o')
    plt.plot(x, 1.0/x, label = '1/N_eff')
    plt.ylabel('expectation value of variance(distance modulus)')
    plt.xlabel('N_effective')
    plt.ylim(-1, 2)
    plt.savefig('preliminary_bump_s_1_3.pdf')
    
    
