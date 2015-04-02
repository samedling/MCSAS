#!/usr/bin/python

import os, sys, pickle, random, pylab, time
import numpy as np
import Image
from scipy.optimize import leastsq

#Contains subroutines to be used by newgui for fitting.

#New User Variables:
#crop=(xmin,xmax,ymin,ymax)
#maxiter (20-whatever; default 0 is rather large)

#1.   Load dictionary from defaults (for now) or gui (ultimately).  New user_variables: crop_x, crop_y (separate min & max?).
#2.   Load TIFF, downsample, apply mask.
#3.   Fit parameters: r1,r2,l,rho1,rho2,rot_z



def load_exp_image(filename,downsample=(100,100),crop=(0,0,0,0),mask_threshold=100):
   '''Loads experimental data from file, cropping it and downsampling it if neccessary, and normalizes.  Also outputs the mask corresponding to the beamstop.'''
   ##If this takes a long time to run, make separate routine to load_exp_data already cropped, downsampled, and in a numpy array.
   #exp_data=np.asarray(Image.open(filename))
   #global downsample,crop
   if sum(crop):
      img=Image.open(filename)
      crop_to=(crop[0],crop[1],im.size[0]-crop[3],im_size[1]-crop[4])      #(left,top,right,bottom), in the image view
      cropped=img.crop(crop_to)
   else:
      cropped=Image.open(filename)
   downsampled=cropped.resize(downsample,Image.ANTIALIAS)      #NEAREST,BILINEAR,BICUBIC,ANTIALIAS (worst to best; fastest to slowest)
   exp_data=numpy.array(downsampled)
   mask=np.ones(product(exp_data.shape)).reshape(exp_data.shape)
   for i in range(mask.shape[0]):
      for j in range(mask.shape[1]):
         if exp_data[i,j] < mask_threshold:        #do after normalize?
            mask[i,j] = 0
   normalize = 1.0/np.sum(exp_data)
   exp_data=exp_data*normalize
   #img=Image.fromarray(exp_data)   #To go back to an image.
   return exp_data,mask

def sync_dict(parameters):
   '''Copies parameters to dictionary_SI.'''
   global dictionary_SI
   parameters = r1,r2,z_dim,rho1,rho2    #z_theta,amplitude?
   dictionary_SI['radius_1'] = r1
   dictionary_SI['radius_2'] = r2
   dictionary_SI['rho_1'] = rho1
   dictionary_SI['rho_2'] = rho2
   dictionary_SI['z_dim'] = z_dim

def residuals(param,exp_data,mask=1):
   '''Returns residual array of difference between experimental data and data calculated from passed parameters.'''
   #global dictionary_SI
   sync_dict(param)
   err = np.zeros(product(exp_data.shape)).reshape(exp_data.shape)
   load_functions()    #DO I NEED?  #Reintilizes functions with the new parameters.
   #calc_intensity = Average_Intensity() #might just take longer or might be necessary to accomodate randomness in Points_For_Calculation
   calc_intensity = Detector_Intensity(Points_For_Calculation())  #like Average_Intensity() but just runs once and without time printouts
   normalize = 1.0/np.sum(calc_intensity)
   for i in range(len(n_pixels)):
      for j in range(len(n_pixels)):
         err[i,j] = exp_data[i,j]-calc_intensity[i,j]*normalize
   return err.ravel(mask*err)  #flattens err since leastsq only takes a 1D array

def fit_data(filename,maxiter=20):
   '''Loads experimental data from filename, fits the data using current dictionary as initial guesses, leaves final parameters in dictionary.'''
   exp_data,mask=load_exp_image(filename)
   #Get initial guess:
   guess=[dictionary_SI['radius_1'],dictionary_SI['radius_2'],dictionary_SI['rho_1'],dictionary_SI['rho_2'],dictionary_SI['z_dim'],]
   fit_param = leastsq(residuals,guess,args=(exp_data),full_output=1,maxfev=maxiter)
   print('Fit converged after {0} function calls.'.format(fit_param[2]['nfev']))
   print(fit_param[0])
   sync_dict(fit_param[0])    #Save Final Fit Parameters
   #plot difference plot or (experimental data and final plot)
   fit=Average_Intensity()
   save(fit,_fit)
   diff=residuals(fit_param[0],exp_data).reshape(exp_data.shape)
   save(diff,_fit_residuals)
   view_fit(fit,diff)

def view_fit(exp_data,fit,residuals):
   '''Copied from view_intensity() with minor changes to plot all relevant fitting plots.'''
   global dictionary,dictionary_SI,dictionary_in,sim_info
   get_numbers_from_gui()
   load_functions()
   Intensity_plot(exp_data,"exp_data",dictionary_SI['title'],0)
   Intensity_plot(fit,"fit",dictionary_SI['title'],1)
   Intensity_plot(residuals,"residuals",dictionary_SI['title'],2)
   clear_mem()
   print "Program Finished"



