#!/usr/bin/python
#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

import os, sys, pickle, random, pylab, time
import numpy as np
from PIL import Image
#from scipy import misc  #alternative to PIL/Image in case that doesn't work.
from scipy.optimize import leastsq

#Contains subroutines to be used by newgui for fitting.
#Fit parameters: r1,r2,z_dim(length),rho1,rho2,x_theta,y_theta,z_theta

#To Do:
#Add ability to import experimental data from files already cropped/downsampled.

#global load_functions,get_numbers_from_gui,Intensity_plot
#from newgui import load_functions,get_numbers_from_gui
#from Plotting_Functions import Intensity_plot

class Fit_Parameters():
   def __init__(self):
      global dictionary,dictionary_SI
      shape=dictionary['shape']
      always=('x_theta','y_theta','z_theta')
      self.names=[]
      if shape in (1,2,6):    #Sphere, Cylinder, or Hex Prism
         self.density_params=('radius_1','rho_1')
      elif shape == 3:        #Core Shell
         self.density_params=('radius_1','radius_2','rho_1','rho_2')
      elif shape == 4:        #Gaussian
         self.density_params=('radius_2')
      elif shape in (5,14):   #Chopped Cone, Double Done
         self.density_params=('radius_1','radius_2','rho_1','z_dim')
      elif shape == 7:        #Rect. Prism#wHY NOT WITH 1,2,6??
         self.density_params=('radius_2','rho_1')
      elif shape in (8,11):   #Bubbles, Double Slit
         self.density_params=('radius_1','radius_2','rho_1')
      elif shape in (9,12,13):#Chopped Cylinder, N-Shaped Chopped Cone, Sine
         self.density_params=('radius_1','radius_2','rho_1','rho_2','z_dim')
      elif shape == 10:
         print('Model not supported.')
      else:
         print('Unknown model. Assuming model uses all parameters.')
         self.density_params=('radius_1','radius_2','rho_1','rho_2','z_dim')
      for name in (self.density_params+always):
         if dictionary['fit_'+name]:
            self.names.append(name)
      self.values=[dictionary_SI[var] for var in self.names]
      self.length=len(self.values)
      self.units=[]
      for i in range(self.length):
         if self.names[i][2:] == 'theta':
            if dictionary['degrees']:
               self.units.append(' degrees')
            else:
               self.units.append(' radians')
         elif self.names[i][:-2] == 'radius':
            self.units.append(' nm')
         elif self.names[i] == 'z_dim':
            self.units.append(' nm')
         else:
            self.units.append('')

   def sync_dict(self):
      '''Copies parameters to dictionary_SI.'''
      global dictionary_SI
      for i in range(self.length):
         dictionary_SI[self.names[i]] = self.values[i]

   def print_param(self):
      global dictionary
      convert_from_SI()
      for i in range(self.length):
         print('{0} is {1}{2}.'.format(self.names[i],self.values[i],self.units[i]))

   def get_param(self):
      return self.values

   def set_param(self,parameters):
      self.values = parameters

   #def update_gui(self):




#def load_exp_image():

def convert_from_SI():
    global dictionary,dictionary_SI
    dictionary["radius_1"] = dictionary_SI["radius_1"]*10**9
    dictionary["radius_2"] = dictionary_SI["radius_2"]*10**9
    dictionary["z_dim"] = dictionary_SI["z_dim"]*10**9
    if dictionary["degrees"] == 1: #Conveting from radians
       dictionary["x_theta"] = dictionary_SI["x_theta"]*180/np.pi
       dictionary["y_theta"] = dictionary_SI["y_theta"]*180/np.pi
       dictionary["z_theta"] = dictionary_SI["z_theta"]*180/np.pi


