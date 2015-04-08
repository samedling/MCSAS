#!/usr/bin/python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import os, sys, pickle, random, pylab, time
import numpy as np
from PIL import Image
#from scipy import misc  #alternative to PIL/Image in case that doesn't work.
from scipy.optimize import leastsq

#Contains subroutines to be used by newgui for fitting.
#Fit parameters: r1,r2,z_dim(length),rho1,rho2,z_theta

#To Do:
#Figure out image importing on Mac OS X?
#Convert SI parameters back to non-SI.
#Need to run make SI thing?
#Add ability to import experimental data from files already cropped/downsampled.
#Allow to crop by specifying center instead of edges?
#Add checkboxes to hold some fit parameters fixed.



def load_exp_image():
   '''Loads experimental data from file, cropping it and downsampling it if neccessary, and normalizes.  Also outputs the mask corresponding to the beamstop.'''
   global dictionary
   downsample=(dictionary['pixels'],dictionary['pixels'])
   crop=dictionary['crop']
   mask_threshold=dictionary['mask_threshold']
   filename=dictonary['filename']
   #exp_data=np.asarray(Image.open(filename))
   if sum(crop):
      img=Image.open(filename)
      crop_to=(crop[0],crop[1],img.size[0]-crop[3],img_size[1]-crop[4])      #(left,top,right,bottom), in the image view
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

def plot_exp_data():
    global dictionary_SI
    image=load_exp_image()
    Intensity_plot(image[0],"exp_data",dictionary_SI['title'],0)
    return

def sync_dict(parameters):
   '''Copies parameters to dictionary_SI.'''
   global dictionary_SI
   parameters = r1,r2,z_dim,rho1,rho2,z_theta
   dictionary_SI['radius_1'] = r1
   dictionary_SI['radius_2'] = r2
   dictionary_SI['z_dim'] = z_dim
   dictionary_SI['rho_1'] = rho1
   dictionary_SI['rho_2'] = rho2
   dictionary_SI['z_theta'] = z_theta

def print_parameters():
   global dictionary_SI
   print('Radius 1 is {0}.'.format(dictionary_SI['radius_1']))
   print('Radius 2 is {0}.'.format(dictionary_SI['radius_2']))
   print('Length is {0}.'.format(dictionary_SI['z_dim']))
   print('Rho 1 is {0}.'.format(dictionary_SI['rho_1']))
   print('Rho 2 is {0}.'.format(dictionary_SI['rho_2']))
   print('z rotation is {0}.'.format(dictionary_SI['z_theta']))
   return

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

def fit_step(exp_data,update_freq=20):
   '''Runs a small number of iterations of fitting exp_data.'''
   global dictionary_SI
   guess=[dictionary_SI['radius_1'],dictionary_SI['radius_2'],dictionary_SI['z_dim'],dictionary_SI['rho_1'],dictionary_SI['rho_2'],dictionary['z_theta']]
   fit_param = leastsq(residuals,guess,args=(exp_data),full_output=1,maxfev=update_freq)
   sync_dict(fit_param[0])    #Save Final Fit Parameters
   return fit_param

def perform_fit():
   '''Loads experimental data from filename, fits the data using current dictionary as initial guesses, leaves final parameters in dictionary.'''
   global dictionary,dictionary_SI
   filename = dictionary['fit_file']
   maxiter = dictionary['max_iter']
   update_freq = dictionary['update_freq']
   #need to refresh dictionary_SI?
   total_steps = 0
   if update_freq == 0:
      update_freq = max_iter
   exp_data,mask=load_exp_image()
   print('{0}: Starting fit...'.format(time.strftime("%X")))
   total_steps = 0
   while total_steps < maxiter:
      fit_param = fit_step(exp_data,update_freq)
      total_steps+=fit_param[2]['nfev']
      if fit_param[2]['nfev'] < update_freq:      #Better parameter to see if fit is completed?
         print('{0}: Converged after {1} function calls.'.format(time.strftime("%X"),total_steps))
         print_parameters()
         break
      else: #This line not necessary, but I think it improves readability.
         print('{0}: On function call {1}...'.format(time.strftime("%X"),total_steps))
         print(fit_param[0])
         view_fit(exp_data,fit,diff)
   if total_steps >= maxiter:
      print('{0}: Fit did not converge in {1} steps.'.format(time.strftime("%X"),total_steps))
   fit=Average_Intensity()
   save(fit,_fit)
   #diff=residuals(fit_param[0],exp_data).reshape(exp_data.shape)
   diff=fit_param[2]['fvec'].reshape(exp_data.shape)
   save(diff,_fit_residuals)
   view_fit(exp_data,fit,diff)

def view_fit(exp_data,fit,residuals):
   '''Copied from view_intensity() with minor changes to plot all relevant fitting plots.'''
   global dictionary,dictionary_SI
   plot_fit = dictionary['plot_fit_tick']
   plot_residuals = dictionary['plot_residuals_tick']
   #get_numbers_from_gui()
   #load_functions()
   if plot_fit:
      Intensity_plot(exp_data,"exp_data",dictionary_SI['title'],0)
      Intensity_plot(fit,"fit",dictionary_SI['title'],1)
   if plot_residuals:
      Intensity_plot(residuals,"residuals",dictionary_SI['title'],2)
   clear_mem()
   print "Program Finished"





#def fit_data(filename,maxiter=20):
#   '''Loads experimental data from filename, fits the data using current dictionary as initial guesses, leaves final parameters in dictionary.'''
#   exp_data,mask=load_exp_image(filename)
#   print('Starting fit at {0}.'.format(time.strftime("%X")))
#   #Get initial guess:
#   guess=[dictionary_SI['radius_1'],dictionary_SI['radius_2'],dictionary_SI['rho_1'],dictionary_SI['rho_2'],dictionary_SI['z_dim']]
#   fit_param = leastsq(residuals,guess,args=(exp_data),full_output=1,maxfev=maxiter)
#   print('Fit converged at {0} after {1} function calls.'.format(time.strftime("%X"),fit_param[2]['nfev']))
#   print(fit_param[0])
#   sync_dict(fit_param[0])    #Save Final Fit Parameters
#   #plot difference plot or (experimental data and final plot)
#   fit=Average_Intensity()
#   save(fit,_fit)
#   #diff=residuals(fit_param[0],exp_data).reshape(exp_data.shape)
#   diff=fit_param[2]['fvec']
#   save(diff,_fit_residuals)
#   view_fit(exp_data,fit,diff)

#Scipy implementation instead of PIL/Image. Unfinished downsampling.  Unnecessary?
#def load_exp_image_scipy(filename,downsample=(100,100),crop=(0,0,0,0),mask_threshold=100):
#   '''Loads experimental data from file, cropping it and downsampling it if neccessary, and normalizes.  Also outputs the mask corresponding to the beamstop.'''
#   ##If this takes a long time to run, make separate routine to load_exp_data already cropped, downsampled, and in a numpy array.
#   #exp_data=np.asarray(Image.open(filename))
#   #global downsample,crop
#   if sum(crop):
#      img = misc.imread(filename)
#      cropped = img[crop[0]:-crop[2],crop[1]:-crop[3]]
#   else:
#      cropped = misc.imread(filename)
#   downsampled =  #downsample using scipy.interpolate
#   mask=np.ones(product(exp_data.shape)).reshape(exp_data.shape)
#   for i in range(mask.shape[0]):
#      for j in range(mask.shape[1]):
#         if exp_data[i,j] < mask_threshold:        #do after normalize?
#            mask[i,j] = 0
#   normalize = 1.0/np.sum(exp_data)
#   exp_data=exp_data*normalize
#   #img=Image.fromarray(exp_data)   #To go back to an image.
