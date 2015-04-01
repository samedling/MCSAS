#!/usr/bin/python

import os, sys, pickle, random, pylab, time
import numpy as np
from scipy.optimize import leastsq

#Contains subroutines to be used by newgui for fitting.

#New User Variables:
#crop=(xmin,xmax,ymin,ymax)
#maxiter (20-whatever; default 0 is rather large)

#1.   Load dictionary from defaults (for now) or gui (ultimately).  New user_variables: crop_x, crop_y (separate min & max?).
#2.   Load TIFF, downsample, apply mask.
#3.   Fit parameters: r1,r2,l,rho1,rho2,rot_z



def load_exp_data(filename,crop=(0,0,0,0),downsample=(100,100)):     #UNFINISHED
   '''Loads experimental data from file, cropping it and downsampling it if neccessary, and normalizes.'''
   #downsampling might not need argument; can use global dictionary
   #Also return mask?
   return exp_data

def sync_dict(parameters):
   '''Copies parameters to dictionary_SI.'''
   global dictionary_SI
   parameters = r1,r2,z_dim,rho1,rho2    #z_theta,amplitude?
   dictionary_SI['radius_1'] = r1
   dictionary_SI['radius_2'] = r2
   dictionary_SI['rho_1'] = rho1
   dictionary_SI['rho_2'] = rho2
   dictionary_SI['z_dim'] = z_dim

def residuals(param,exp_data):
    #global dictionary_SI
    sync_dict(param)
    err = np.zeros(product(exp_data.shape)).reshape(exp_data.shape)
    #load_functions()    #DO I NEED?  #Reintilizes functions with the new parameters.
    calc_intensity = Detector_Intensity(Points_For_Calculation())  #calc_intensity = Average_Intensity() might just take longer or might be necessary to accomodate randomness in Points_For_Calculation
    normalize = 1.0/np.sum(calc_intensity)
    for i in range(len(n_pixels)):
        for j in range(len(n_pixels)):
            err[i,j] = exp_data[i,j]-calc_intensity[i,j]*normalize      #Secaling factor and/or mask?
    return err.ravel(err)  #flattens err since leastsq only takes a 1D array

def fit_data(filename,maxiter=20):
   '''Loads experimental data from filename, fits the data using current dictionary as initial guesses, leaves final parameters in dictionary.'''
   exp_data=load_exp_data(filename)
   #Get initial guess:
   guess=[dictionary_SI['radius_1'],dictionary_SI['radius_2'],dictionary_SI['rho_1'],dictionary_SI['rho_2'],dictionary_SI['z_dim'],]
   fit_param = leastsq(residuals,guess,args=(exp_data),full_output=1,maxfev=maxiter)
   print('Fit converged after {0} function calls.'.format(fit_param[2]['nfev']))
   #print(fit_param[0])
   sync_dict(fit_param[0])    #Save Final Fit Parameters
   #plot difference plot or (experimental data and final plot)
   fit=Average_Intensity()
   save(fit,_fit)
   diff=residuals(fit_param[0],exp_data).reshape(exp_data.shape)
   save(diff,_fit_residuals)
   #make viewintensity show things, or even better, make a viewfit to show everything at once!


