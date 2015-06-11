#!/usr/bin/python
version = '0.3.2'


try:
   from Tkinter import *
except:
   from tkinter import *

import os, sys, pickle, random, pylab, time
import numpy as np
from scipy.special import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

from PIL import Image
from scipy.optimize import leastsq
#from scipy import misc     #alternative to PIL
#from scipy import ndimage  #possible smoothing of exp_data before viewing

import global_vars as g

#Looks for fastmath.so to speed up intensity calculation.
if g.opencl_enabled:
   try:
      import pyopencl as cl
      from opencl import OpenCL
      g.opencl_density = OpenCL()
      g.opencl_density.load_program('density.cl')
      g.opencl_sumint = OpenCL()
      g.opencl_sumint.load_program('sumint.cl')
      print("Accelerating using OpenCL.")
   except ImportError:
      g.opencl_enabled = False
if g.f2py_enabled:
   try:
      import fastmath
      print("Accelerating using f2py.")
   except ImportError:
      g.f2py_enabled = False
if not g.f2py_enabled and not g.opencl_enabled:
   print("Could not accelerate using either OpenCL or f2py.")
   print("See README for how to install either OpenCL or f2py.")
   print("In the meantime, fitting is not recommended.")

if g.debug:
   np.random.seed([2015])     #Locks random seed to allow for speedtesting.
   g.verbose = True


#These are the default settings
g.dictionary = {'advanced':1, 'altitude':45, 'analytic': 2, 'ave_dist': 0.6, 'azimuth':45, 'bound': 1, 'circ_delta':5, 'comments':'',
              'degrees': 1, 'energy_wavelength': 12, 'energy_wavelength_box': 0, 'gauss':0, 'log_scale': 1, 'maximum': 0.01, 'minimum': 1e-8,
              'num_plots': 1, 'pixels': (200,200), 'proportional_radius':0.5, 'QSize': 6,'Qz': 0, 'radius_1': 5.0, 'radius_2': 2.5, 'rho_1': 1.0, 'rho_2': -0.5,
              'save_img':1, 'save_name': 'save_name', 'scale': 1,'SD':1, 'seq_hide':0, 'shape': 2, 's_start': 0, 's_step': 2,
              's_stop': 1, 'subfolder':'subfolder', 's_var': 'x_theta', 'symmetric': 0,
              'theta_delta':20, 'ThreeD': 0, 'title': 'title', 'x_theta': 0,'y_theta': 0,'z_theta': 0,'z_dim': 100,'z_scale':1,#}
              'fit_file': 'fit_file', 'center': (0,0), 'border': 0, 'max_iter': 1000, 'update_freq': 0, 'plot_fit_tick': 1, 'plot_residuals_tick': 1, 'mask_threshold': 10, 'background': 2e-5, 'grid_compression': 5,
              'fit_radius_1': 1, 'fit_radius_2': 0, 'fit_rho_1': 1, 'fit_rho_2': 0, 'fit_z_dim': 1, 'fit_x_theta': 1, 'fit_y_theta': 1, 'fit_z_theta': 1, 'fit_background': 1, 'fit_other': 0
              }

#####            Importing data or using defaults              #############

if os.name == 'nt':  #or sys.platform == 'win32'
   root_folder = os.path.dirname(sys.argv[0]) #Windows
else: #os.name == 'posix' or sys.platform == 'linux2'
   root_folder = os.getcwd()  #Mac/Linux
#Check for write access?

#Windows: os.name='nt', sys.platform='win32'
#Ubuntu: os.name='posix', sys.platform='linux2'
#Mac OS X: os.name='posix', sys.platform='darwin'

try:
    d = pickle.load(open(root_folder+"/default.txt", 'rb'))
    if len(g.dictionary) != len(d): #I check that it is the same length - This is needed if any new variables are added to g.dictionary
        a= 1/0
    g.dictionary = d
except:
    print "Previously used variables could not be loaded. \nUsing default settings instead."
    with open(root_folder+"/default.txt", 'wb') as f:
       pickle.dump(g.dictionary, f)

g.dictionary = {x:g.dictionary[x] for x in g.dictionary} #This contains the unaltered parameters
g.dictionary_in = {a:g.dictionary[x] for x in g.dictionary for a in [x, x+'2']} #This contains raw data from the GUI. The +'2' is so that i can control checkboxes
g.dictionary_SI = {x:g.dictionary[x] for x in g.dictionary} #this g.dictionary has the parameters after they have been converted to SI units

#This is the list of all the Monte Carlo Models that you choose from.
MC_num_and_name = np.array([["Analytic Model Only",0],
                        ["Sphere",1],
                        ["Cylinder",2],
                        ["Core shell cylinder",3],
                        ["Gaussian",4],
                        ["Cone Model",5],
                        ["Hexagonal Prism",6],
                        ["Rectangular Prism",7],
                        ["String of Bubbles",8],
                        ["Chopped up Cylinder",9],
                        ["Custom CSV Defining the Radius",10],
                        ["Double Slit",11],
                        ["N-gon Truncated Cone",12],
                        ["Sine Shaped Oscillation",13],
                        ["Double Cone",14],
                        ["Eliptical Cylinder",15]
                        ])
MC_num_and_name_dict = {x[0]:x[1] for x in MC_num_and_name} #This is needed, so that when an option is chosen, we can find the shape number.

#This is the list of all the analytic models that you choose from.
Analytic_options = np.array([["None",0],
                        ["Sphere",1],
                        ["Cylinder",2],
                        ["Core shell cylinder",3],
                        ["Gaussian - check this formula",4]
                        ])
Analytic_dict = {x[0]:x[1] for x in Analytic_options} #This is needed, so that when an option is chosen, we can find the shape number.


def xy_dim(): #Since x_dim and y_dim are not defined by the user any more, the function to define them is here.
   g.dictionary_SI['x_dim'] = 2.*g.dictionary_SI['radius_1']
   g.dictionary_SI['y_dim'] = 2.*g.dictionary_SI['radius_1']
   
   

def get_numbers_from_gui():
    #Here I get all parameters from the GUI and put them into g.dictionary
    for x in g.dictionary:
       if x=='comments':
          g.dictionary[x] = g.dictionary_in[x].get(1.0,END).rstrip()
       elif x=='advanced' or x== 'seq_hide':
          g.dictionary[x] = g.dictionary[x]
       elif x=='shape':
          g.dictionary[x] = MC_num_and_name_dict[g.dictionary_in['shape'].get()]
       elif x=='analytic':
          g.dictionary[x] = Analytic_dict[g.dictionary_in['analytic'].get()]
       else:
          try:
             g.dictionary[x] = g.dictionary_in[x].get() 
          except AttributeError:
             if g.debug:
                print("{0} could not be imported from GUI.".format(x))
    
    for x in g.dictionary:
        try:
            g.dictionary[x] = float(g.dictionary[x]) #I am turning the numbers from strings to floats.
            if g.dictionary[x]==int(g.dictionary[x]):#I am turning floats to integers if they are the same
               g.dictionary[x] = int(g.dictionary[x])
        except:
            None
    if not os.path.exists(root_folder+'/'+g.dictionary_SI['subfolder']):#making the subfolder, if it doesn't exist
       os.makedirs(root_folder+'/'+g.dictionary_SI['subfolder'])
       time.sleep(2) #Making the subfolder takes a few seconds, so we need to delay the program, otherwise it will try save things into the folder before it is made.

    g.dictionary_SI = {x: g.dictionary[x] for x in g.dictionary}
    g.dictionary_SI['path_to_subfolder'] = os.path.join(root_folder,g.dictionary['subfolder'],g.dictionary['save_name']) #This is for convienience

    #Converting to SI units.
    g.dictionary_SI["z_dim"] = g.dictionary["z_dim"]*10**-9
    g.dictionary_SI["ave_dist"] = g.dictionary["ave_dist"]*10**-9
    g.dictionary_SI["travel"] = g.dictionary_SI["ave_dist"]
    g.dictionary_SI["radius_1"] = g.dictionary["radius_1"]*10**-9
    xy_dim()#defining x_dim and y_dim - dependent of radius_1
    g.dictionary_SI["radius_2"] = g.dictionary["radius_2"]*10**-9
    g.dictionary_SI["QSize"] = g.dictionary["QSize"]*10**9

    #g.dictionary_SI['num_plot_points'] = int(g.dictionary_SI['pixels']/2.)
    g.dictionary_SI['num_plot_points'] = min([int(i) for i in g.dictionary['pixels'].split()])/2
    g.dictionary_SI['delta'] = 1. #number of pixels in width

    g.dictionary_SI['circ_delta'] = 1.4*g.dictionary_SI['circ_delta']
    try:
       g.dictionary_SI['pixel_radius'] = g.dictionary['proportional_radius']*g.dictionary_SI['pixels']/2.
    except TypeError:
       g.dictionary_SI['pixel_radius'] = g.dictionary['proportional_radius']*int(g.dictionary_SI['pixels'].split()[0])/2.
    g.dictionary_SI['theta_delta'] = 6.283/g.dictionary['theta_delta']
    


    if g.dictionary_SI["energy_wavelength_box"] == 0: #Checkbox
       g.dictionary_SI["EHC"] = 2.*np.pi*(g.dictionary_SI["energy_wavelength"])*1.602176487*10**10/(6.62606896*2.99792458) #Energy to 2pi/lambda
    else:
       g.dictionary_SI["EHC"] = 2. * np.pi / g.dictionary_SI["energy_wavelength"]#lambda to 2pi/lambda
    if g.dictionary_SI["degrees"] == 1: #Conveting to radians
       g.dictionary_SI["x_theta"] = g.dictionary["x_theta"]*np.pi/180
       g.dictionary_SI["y_theta"] = g.dictionary["y_theta"]*np.pi/180
       g.dictionary_SI["z_theta"] = g.dictionary["z_theta"]*np.pi/180

    with open(root_folder+"/default.txt", 'wb') as f:
        pickle.dump(g.dictionary, f)#Saving the infomation from g.dictionary so it can be loaded later

    with open(g.dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
       pickle.dump(g.dictionary, f) #saving a copy in the subfolder for reference.


def load_functions(): #This loads the functions from the other files. It needs to be dynamic, hence I cannot use import.
   execfile(root_folder+"/Monte_Carlo_Functions.py",globals())
   execfile(root_folder+"/Plotting_Functions.py", globals())
   execfile(root_folder+"/density_formula.py", globals())
   execfile(root_folder+"/analytic_formula.py", globals())

def change_units(number): #Used for sequences. A value is converted to SI units.
        for x in g.dictionary_SI:
           if x == g.dictionary_SI['s_var']:
              g.dictionary_SI[x] = number*10**-9
              if x == 'ave_dist':
                 g.dictionary_SI['travel'] = g.dictionary_SI[x]
              if x == 'travel':
                 g.dictionary_SI['ave_dist'] = g.dictionary_SI[x]
              if g.dictionary_SI['s_var'] != 'y_dim' and g.dictionary_SI['s_var'] != 'x_dim':
                 xy_dim()#This is here, mainly for the double slit - if you want to make the slits higher, you can. Most other functions are radially symmetric.
                               
              if x == "QSize":
                 g.dictionary_SI[x] = number*10**9
              if x == "energy_wavelength":
                  if g.dictionary_SI['energy_wavelength_box'] == 1:
                     g.dictionary_SI['EHC'] = 2.*np.pi*number*1.602176487*10**10/(6.62606896*2.99792458)
                  else:
                     g.dictionary['EHC'] = 2. * np.pi / number
              if x == "x_theta" or x == "y_theta" or x == "z_theta":
                  if g.dictionary_SI['degrees'] == 1:
                      g.dictionary_SI[x] = number*np.pi/180

def save_vars_to_file(extra): #here I save all the infomation into a text file that is easy to read. extra is a string of extra infomation that you might want to include.
   get_numbers_from_gui()
   sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
   sim_info.write("\n\n\nDATE: "+time.strftime("%d")+time.strftime("%b")+time.strftime("%y")+"    TIME: "+time.strftime("%X")+"\n")
   sim_info.write(g.dictionary_SI['comments'])
   print g.dictionary_SI['path_to_subfolder']
   print g.dictionary['comments']
   print extra
   sim_info.write(extra+"\n")
   for x in sorted(g.dictionary, key=lambda x: x.lower()):
      if x != 'comments':
         sim_info.write(x+" : " +str(g.dictionary[x])+'\n')
   sim_info.close()
   
def clear_mem():#This function clears memory to try and reduce mem useage.
   try:                     #You don't actually need try/except, but if you are,
      Intensity = None      #then you should "except NameError", not everything.
   except:
      None
   try:
      points = None
   except:
      None
   try:
      radial_intensity = None
   except:
      None
   try:
      cumulative = None
   except:
      None

###############           Buttons used in the GUI         ##############

#This is a list of some Common Variables for use in sequence.
ALLVARIABLES = [['altitude', 'Altitude'],
                ["ave_dist","Neighbouring Point Distance"],
                ['azimuth', 'Azimuth'],
                ["energy_wavelength","Energy/Wavelength"],
                ["pixels","pixels on the detector"],
                ["QSize","Q Range of the entire detector"],
                ["radius_1","Radius 1"],
                ["radius_2","Radius 2"],
                ["rho_1","Rho 1"],
                ["rho_2","Rho 2"],
                ["x_theta","rotation in x direction"],
                ["y_theta","rotation in y direction"],
                ["z_dim","Length"],
                ['z_scale', 'z-direction scaling of\nneighbouring point distance'],
                ["z_theta","rotation in z direction"],
                ]
def show_sequence_variables(): #Common Variables button, displays ALLVARIABLES, above.
   dens_options = Tk()
   dens_options.title("Sequence options")
   Label(dens_options, text = "Choose a Variable to Change", font = "Times 14 bold underline").grid(row = 0, column = 0, columnspan = 2, sticky = W)
   Label(dens_options, text = "Variable", font = "Times 11 underline").grid(row = 1, column = 0, sticky = W)
   Label(dens_options, text = "Description", font = "Times 11 underline").grid(row = 1, column = 1, sticky = W)
   for ROW in range(len(ALLVARIABLES)):
      description = ALLVARIABLES[ROW]
      Label(dens_options, text = description[0]).grid(row = ROW+3, column = 0, sticky = W)
      Label(dens_options, text = description[1]).grid(row = ROW+3, column = 1, sticky = W)
   dens_options.mainloop()

def plot_points(): #This runs the Real Space to plot the points in Real Space
    get_numbers_from_gui()
    save_vars_to_file("Plot Points")
    load_functions()
    if g.dictionary['seq_hide'] == 1:
       if g.dictionary['gauss']==0:
          current_value = g.dictionary['s_start']
       else:
          current_value = np.random.normal(loc = g.dictionary[g.dictionary['s_var']], scale = g.dictionary['SD'])
       change_units(current_value)
    Points_Plot(Points_For_Calculation(), 'points', 1)
    clear_mem()
    print "Program Finished"

def view_intensity(): #This allows you to view a premade intensity
    get_numbers_from_gui()
    load_functions()
    radial_intensity = pylab.loadtxt(g.dictionary_SI['path_to_subfolder']+"radial_intensity.csv", delimiter=",")
    radial_intensity_plot(radial_intensity, "radial", g.dictionary_SI['title'], 0)
    Intensity = pylab.loadtxt(g.dictionary_SI['path_to_subfolder']+"intensity.csv", delimiter=",")
    Intensity_plot(Intensity, "intensity", g.dictionary_SI['title'], 1)
    clear_mem()
    print "Program Finished"
    
def make_intensity(): #This makes an intensity
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file("Monte Carlo Intensity")
    load_functions()
    Intensity = Average_Intensity()
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    if g.dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    print "Program Finished"
    
def sequence(): #This makes a sequence of intensities
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file("Monte Carlo Sequence")
    load_functions()
    for frame_num in range(int(g.dictionary['s_step'])):
       sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
       sim_info.write("\nFrame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step'])))
       sim_info.close()

       print "\nmaking frame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step']))
       if g.dictionary['gauss']==0:
          try:
             current_value = (g.dictionary['s_stop']-g.dictionary['s_start'])*frame_num/(1.*g.dictionary['s_step']-1.)+1.*g.dictionary['s_start']
          except ZeroDivisionError:
             current_value = g.dictionary['s_start']
       else:
          current_value = np.random.normal(loc = g.dictionary[g.dictionary['s_var']], scale = g.dictionary['SD'])
       change_units(current_value)
       Intensity = Average_Intensity()
       save(Intensity, "intensity"+str(frame_num+1))
       radial_intensity = radial(Intensity)
       save(radial_intensity, "radial_intensity"+str(frame_num+1))
       try:
           cumulative += np.asarray(Intensity)
       except NameError:
           cumulative = np.asarray(Intensity)
       title = g.dictionary_SI['title']+" "+g.dictionary_SI['s_var']+'='+str(current_value)
       if g.dictionary_SI['save_img'] == 1:
          Intensity_plot(Intensity, "intensity" + str(frame_num+1), title, 0)
          radial_intensity_plot(radial_intensity, "radial_intensity" + str(frame_num+1), title, 0)
       clear_mem()
    Intensity = cumulative / g.dictionary_SI['s_step']
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    g.dictionary_SI['title'] = g.dictionary_SI['title']+" Averaged" + g.dictionary_SI['s_var']
    if g.dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    print "Program Finished"


def theory_plot(): #This plots an analytic model
   get_numbers_from_gui()
   save_vars_to_file("Analytic Intensity")
   load_functions()
   Intensity = theory_csv()
   save(Intensity, "intensity")
   radial_intensity = radial(Intensity)
   save(radial_intensity, "radial_intensity")
   if g.dictionary_SI['save_img'] == 1:
      view_intensity()
   clear_mem()
   print "Program Finished"



def theory_seq(): #This plots a sequence created with the analytic model
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file("Analytic Sequence")
    load_functions()
    for frame_num in range(int(g.dictionary['s_step'])):
       sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
       sim_info.write("\nFrame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step'])))
       sim_info.close()

       print "\nmaking frame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step']))
       if g.dictionary['gauss']==0:
          try:
             current_value = (g.dictionary['s_stop']-g.dictionary['s_start'])*frame_num/(1.*g.dictionary['s_step']-1.)+1.*g.dictionary['s_start']
          except ZeroDivisionError:
             current_value = g.dictionary['s_start']
       else:
          current_value = np.random.normal(loc = g.dictionary[g.dictionary['s_var']], scale = g.dictionary['SD'])
       change_units(current_value)
       Intensity = theory_csv()
       save(Intensity, "intensity"+str(frame_num+1))
       radial_intensity = radial(Intensity)
       save(radial_intensity, "radial_intensity"+str(frame_num+1))
       try:
           cumulative += np.asarray(Intensity)
       except NameError:
           cumulative = np.asarray(Intensity)
       title = g.dictionary_SI['title']+" "+g.dictionary_SI['s_var']+'='+str(current_value)
       if g.dictionary_SI['save_img'] == 1:
          Intensity_plot(Intensity, "intensity" + str(frame_num+1), title, 0)
          radial_intensity_plot(radial_intensity, "radial_intensity" + str(frame_num+1), title, 0)
       clear_mem()
    Intensity = cumulative / g.dictionary_SI['s_step']
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    g.dictionary_SI['title'] = g.dictionary_SI['title']+" Averaged" + g.dictionary_SI['s_var']
    if g.dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    print "Program Finished"



def circ(): #This plots a the angle at a fixed radius
   get_numbers_from_gui()
   load_functions()
   Intensity = np.asarray(pylab.loadtxt(g.dictionary_SI['path_to_subfolder']+"intensity.csv", delimiter=","))
   data = plotting_circle(Intensity)
   radial_intensity_plot(data, "theta"+str(g.dictionary['radius_2']), g.dictionary['title']+" "+str(g.dictionary['radius_2']), 0)
   angle_plot(data, "Angle"+str(g.dictionary['radius_2']), g.dictionary['title']+" "+str(g.dictionary['radius_2']), 1)
   print "finsihed"
   

def int_seq(): #This is the button, it runs a sequence or a single image depending on whether or not you can edit a sequence (For both analytic models and Monte Carlo Models)
   if g.dictionary['seq_hide'] == 0:
      if g.dictionary['shape'] ==0:
         theory_plot()
      else:
         make_intensity()
   else:
      if g.dictionary['shape']==0:
         theory_seq()
      else:
         sequence()







### Fitting Functions ###


class Fit_Parameters():
   '''Class to keep track of which parameters are being varied, based on shape and checkboxes.
      Contains functions to synchronize these parameters with global dictionaries.
      Also contains list of units for user-friendly output.'''
   def __init__(self):
      shape=g.dictionary['shape']
      always=('x_theta','y_theta','z_theta','background','z_dim')
      self.names=[]
      if shape in (1,2,6):    #Sphere, Cylinder, or Hex Prism
         self.density_params=('radius_1','rho_1')
      elif shape == 3:        #Core Shell
         self.density_params=('radius_1','radius_2','rho_1','rho_2')
      elif shape == 4:        #Gaussian
         self.density_params=('radius_1','radius_2')     #radius_1 only for scaling.
      elif shape in (5,14):   #Chopped Cone, Double Done
         self.density_params=('radius_1','radius_2','rho_1') #z_dim intrinsic too
      elif shape == 7:        #Rect. Prism
         self.density_params=('radius_1','radius_2','rho_1')
      elif shape in (8,11,15):   #Bubbles, Double Slit, Eliptical Cylinder
         self.density_params=('radius_1','radius_2','rho_1')
      elif shape in (9,12,13):#Chopped Cylinder, N-Shaped Chopped Cone, Sine
         self.density_params=('radius_1','radius_2','rho_1','rho_2')  #z_dim intrinsic too
      elif shape == 10:
         print('Model not supported.')
      else:
         print('Unknown model. Assuming model uses all parameters.')
         self.density_params=('radius_1','radius_2','rho_1','rho_2')  #z_dim intrinsic too
      for name in (self.density_params+always):
         if g.dictionary['fit_'+name]:   #Looks for checkbox values.
            self.names.append(name)
      self.values=[g.dictionary_SI[var] for var in self.names]
      self.length=len(self.values)
      self.units=[]
      for i in range(self.length):     #Makes array of unit names correlated with values (for printing).
         if self.names[i][2:] == 'theta':
            if g.dictionary['degrees']:
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
      '''Copies parameters to g.dictionary_SI.'''
      for i in range(self.length):
         g.dictionary_SI[self.names[i]] = self.values[i]

   def print_param(self,logfile=0):
      '''Prints parameters to the screen in a user-friendly format.'''
      convert_from_SI()
      for i in range(self.length):
         lprint('{0} is {1:.4}{2}.'.format(self.names[i],float(g.dictionary[self.names[i]]),self.units[i]),logfile)
         #print('{0} is {1}{2}.'.format(self.names[i],self.values[i],self.units[i])) #Units always wrong?

   def get_param(self):
      '''Returns array usable by fitting routines.'''
      return self.values

   def set_param(self,parameters):
      '''Sets parameters from array from fitting routine.'''
      self.values = parameters

def load_exp_image(preview=False,enlarge_mask=1):
   '''Loads experimental data from file, cropping it and downsampling it if neccessary, and normalizes.  Also outputs the mask corresponding to the beamstop.'''
   downsample=[int(i) for i in g.dictionary_SI['pixels'].split()]
   center=[int(i) for i in g.dictionary['center'].split()]
   border=int(g.dictionary['border'])
   mask_threshold=g.dictionary['mask_threshold']
   filename=g.dictionary['fit_file']
   if preview:
      try:
         exp_data=np.array(Image.open(filename))
      except IOError:
         print('File {0} does not exist.'.format(filename))
         return
      return normalize(exp_data)
   else:
      if sum(center):
         img=Image.open(filename)
         #crop_to=(crop[0],crop[1],img.size[0]-crop[2],img.size[1]-crop[3])      #(left,top,right,bottom) in the image view, but (left,bottom,right,top) in the matplotlib plot
         radius=min(center[0],center[1],img.size[0]-center[0],img.size[1]-center[1])
         l=center[0]-radius+border
         b=center[1]-radius+border
         r=center[0]+radius-border
         t=center[1]+radius-border
         crop_to=(l,b,r,t)
         cropped=img.crop(crop_to)
         print("Cropped to {0}.".format(cropped.size))
      else:
         cropped=Image.open(filename)
      downsampled=cropped.resize((max(downsample),max(downsample)),Image.BICUBIC)      #NEAREST,BILINEAR,BICUBIC,ANTIALIAS (worst to best; fastest to slowest; except ANTIALIAS does weird things sometimes)
      if downsample[0] < downsample[1]: #tall rectangle
         w,h = downsample
         d = h-w
         downsampled=downsampled.crop((d/2,0,h-d/2,h))
      elif downsample[1] < downsample[0]:  #wide rectangle
         w,h = downsample
         d = w-h
         downsampled=downsampled.crop((0,d/2,w,w-d/2))
      print("Resized to {0}.".format(downsample))
      exp_data=np.array(downsampled)
      padded=np.lib.pad(exp_data,((1,1),(1,1)),'edge')   #pads the array for enlarging mask
      mask=np.ones_like(exp_data)
      for i in range(mask.shape[0]):
         for j in range(mask.shape[1]):
            if exp_data[i,j] < mask_threshold:        #todo? do after normalize to use shown value instead of raw value?
               mask[i,j] = 0
            elif enlarge_mask and padded[i:i+3,j:j+3].min() < mask_threshold:        #do after normalize?
               mask[i,j] = 0  #could set to ~0.1 if want to decrease but not zero it.
      #exp_data=normalize(exp_data,mask)
      #img=Image.fromarray(exp_data)   #To go back to an image.
      return normalize(exp_data,mask),mask

def normalize(data,mask=[],background=0):
   '''Normalizes, taking mask into account and, if neccessary, adding constant background first.  Be careful not to run this twice in a row, or you'll add double the background.'''
   if not len(mask):
      mask = np.ones_like(data)
   if background:
      norm_to = (1.0 - g.dictionary_SI['background']*np.sum(mask))
      if norm_to < 0:   #THIS IS BAD
          print("Background is too high!")
          norm_to = 1.0
      if g.debug:
         print("Normalizing to {0} so will be normalized to 1 when background of {1} is added.".format(norm_to,g.dictionary_SI['background']))
      #total = np.sum(data*mask) + g.dictionary_SI['background']*np.sum(mask)
      total = norm_to/np.sum(np.maximum(data*mask,np.zeros_like(data))) #maximum is to avoid counting negative points
      #total = norm_to/np.sum(data*mask)
      normalized = np.maximum(data*total*mask,np.zeros_like(data)) + g.dictionary_SI['background']*mask  #maximum zeros out negative points
      #data *= total                               #todo: might be faster
      #data += g.dictionary_SI['background']*mask  #todo: might be faster
   else:
      total = 1.0/np.sum(np.maximum(data*mask,np.zeros_like(data)))
      normalized = np.maximum(data*total*mask,np.zeros_like(data))
      #data *= total    #todo: should be a little faster
   if g.debug: 
      #print("Normalized so total value is {0} and lowest value is {1}.".format(np.sum(normalized),np.min(normalized)))
      print("Normalized so total value is {0} and lowest value is {1}.".format(np.sum(normalized*mask),np.min(normalized)))
   return normalized

def plot_exp_data():#threshold=1e-7,zero_value=1e-7):
    '''Plots experimental data, using center and crop parameters if nonzero.  Also shows the effect of grid compression if enabled.'''
    get_numbers_from_gui()
    load_functions()    #Needed for plotting routines.
    if g.dictionary['center'] == "0 0":
        image=load_exp_image(preview=True)
        print('Original image size is {0} x {1} pixels.'.format(image.shape[0],image.shape[1]))
        print('This takes a minute...')
        threshold=np.median(image)/10     #These 3 lines not necessary since upper/lower bounds is similar
        zero_value=threshold              #but these work automatically fairly well
        image[image<threshold]=zero_value #so they're worth keeping for now.
        Intensity_plot(image,"exp_data",'Before Cropping or Downsampling',1)
    else:
        cropped=load_exp_image()
        masked=cropped[0]*cropped[1]
        #for i in (2,5,10,20,50,80,90,95,98,99):
           #print('{0}th percentile value is {1:0.4}.'.format(i,np.percentile(masked,i)))
        print('Recommended value for background noise parameter is {0:0.4}.'.format(np.percentile(masked,25)))
        if g.dictionary['grid_compression'] > 1:
           fast_mask(cropped[0],cropped[1],g.dictionary['grid_compression'])
           masked=cropped[0]*cropped[1]
           Intensity_plot(masked,"exp_data2",'After Cropping and Downsampling',1)
        else:
           Intensity_plot(masked,"exp_data2",'After Cropping and Downsampling',1)
    return

def residuals(param,exp_data,mask=[],random_seed=2015):
   '''Returns residual array of difference between experimental data and data calculated from passed parameters.'''
   global parameters
   parameters.set_param(param)
   parameters.sync_dict()
   x=range(exp_data.shape[0])
   y=range(exp_data.shape[1])
   if not len(mask):
      mask = np.ones(exp_data.shape)
   calc_intensity = normalize(Detector_Intensity(Points_For_Calculation(seed=random_seed),mask),mask,True)
   #TODO: it shouldn't be able to set the background too high....???
   #TODO: remove +backgrund from next line??
   err = mask*(exp_data - (calc_intensity + g.dictionary_SI['background']))
   #calc_intensity += g.dictionary_SI'background'] - exp_data  #todo: might be faster
   #calc_intensity *= mask                                     #todo: might be faster
   print('{0}: Total error = {1:.4}; sum of squares = {2:.4}'.format(time.strftime("%X"),np.abs(err).sum(),np.square(err).sum()))
   return np.ravel(err)     #flattens err since leastsq only takes 1D array

def lprint(text,filename=0):
   '''Simple routine for logging text printed to the screen.'''
   print(text)
   if filename:
      with open(filename,'a') as log:
         log.write(text+"\n")

def fit_step(exp_data,mask=[],update_freq=50):
   '''Runs a small number of iterations of fitting exp_data.'''
   global parameters
   guess = parameters.get_param()
   #norm_exp_data = exp_data/np.sum(exp_data*mask) #data normalized taking mask into account.
   exp_data = normalize(exp_data,mask)
   fit_param = leastsq(residuals,guess,args=(exp_data,mask,int(time.time())),full_output=1,maxfev=update_freq)
   with open(root_folder+"/default.txt", 'wb') as f:
      pickle.dump(g.dictionary, f)#Saving the infomation from g.dictionary so it can be loaded later
   with open(g.dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
      pickle.dump(g.dictionary, f) #saving a copy in the subfolder for reference.
   return fit_param

def perform_fit():  #Gets run when you press the Button.
   '''Loads experimental data from filename, fits the data using current g.dictionary as initial guesses, leaves final parameters in g.dictionary.'''
   global parameters
   get_numbers_from_gui()
   load_functions()
   filename = g.dictionary['fit_file']
   max_iter = g.dictionary['max_iter']
   update_freq = g.dictionary['update_freq']
   plot_fit=g.dictionary['plot_fit_tick']
   plot_diff=g.dictionary['plot_residuals_tick']
   logfile=g.dictionary_SI['path_to_subfolder']+'fitlog.txt'   #g.dictionary['fitlog']
   grid_compression=g.dictionary['grid_compression']
   if not g.f2py_enabled and not g.opencl_enabled:
      print('Fortran acceleration is NOT enabled!')
      if grid_compression > 1:
         print('Grid compression does not work without Fortran.')
      print('This will probably take a REALLY LONG time.')
   if g.quiet:
      initial_quiet=True
   else:
      initial_quiet=False
   g.quiet = True
   total_steps = 0
   if update_freq == 0:
      update_freq = max_iter
   exp_data,mask=load_exp_image()
   parameters=Fit_Parameters()  #Creates object of parameters with values and names.
   lprint('Starting values:',logfile)
   parameters.print_param(logfile)
   lprint('{0}: Starting fit...'.format(time.strftime("%X")),logfile)
   if grid_compression > 1:
      mask = fast_mask(exp_data,mask,grid_compression)
   total_steps = 0
   while total_steps < max_iter or max_iter == 0:     #TODO: debug or disable update_interval
      fit_param = fit_step(exp_data,mask,update_freq)
      total_steps+=fit_param[2]['nfev']
      if fit_param[2]['nfev'] < update_freq:      #Checks if fit is completed.
         lprint('{0}: Converged after {1} function calls.'.format(time.strftime("%X"),total_steps),logfile)
         lprint('Final parameter values are:',logfile)
         parameters.print_param(logfile)
         break
      elif total_steps < max_iter:
         lprint('{0}: On function call {1}...'.format(time.strftime("%X"),total_steps),logfile)
         lprint('Current parameter values are:',logfile)
         parameters.print_param(logfile)
         #Intensity_plot(fit_param[2]['fvec'].reshape(exp_data.shape),"residuals",'Difference Plot',1)
      else:
         lprint('{0}: Fit did not converge in {1} steps.'.format(time.strftime("%X"),total_steps),logfile)
         lprint('Current parameter values are:',logfile)
         parameters.print_param(logfile)
   diff=fit_param[2]['fvec'].reshape(exp_data.shape)
   g.quiet=initial_quiet
   save(diff,"_fit_residuals")
   #with open(root_folder+"/default.txt", 'wb') as f:
   #   pickle.dump(g.dictionary, f)#Saving the infomation from g.dictionary so it can be loaded later
   #with open(g.dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
   #   pickle.dump(g.dictionary, f) #saving a copy in the subfolder for reference.
   if plot_fit and plot_diff:
      fit_results=normalize(Average_Intensity(),background=True)
      save(fit_results,"_fit")
      diff = abs(exp_data - fit_results)
      #diff = abs(exp_data - (fit_results + g.dictionary_SI['background']))
      Fit_plot(exp_data*mask,fit_results,diff)
   elif plot_fit:
      fit_results=Average_Intensity()
      save(fit_results,"_fit")
      print('Plotting difference.')
      Intensity_plot(fit_results,"residuals",'Difference Plot',1)
   elif plot_diff:
      print('Plotting difference.')
      Intensity_plot(diff,"residuals",'Difference Plot',1)


def fast_mask(exp_data,mask,speedup=5):
   '''Speeds calculation by adding points to mask.  Speedup can be (2,5,10). Returns new mask, but should also edit the old one.'''
   if speedup < 2:
      return mask
   elif speedup == 2:
      mod=2
      percentile=90
      pad=1
   elif speedup == 5:
      mod=3
      percentile=95
      pad=1
   elif speedup == 10:
      mod=5
      percentile=96
      pad=1
   else:
      percentile=min(99,max(50,100-50./speedup))
      mod=speedup
      pad=1
   threshold=np.percentile(exp_data,percentile)
   total=np.product(mask.shape)
   starting=total-mask.sum()
   if pad:
      padded=np.lib.pad(exp_data,((1,1),(1,1)),'edge')   #pads the array
      for i in range(exp_data.shape[0]):
         for j in range(exp_data.shape[1]):
            if padded[i:i+3,j:j+3].max() < threshold:
               if not (i%mod==mod/2 and j%mod==mod/2):
                  mask[i,j] = 0
   else:
      for i in range(exp_data.shape[0]):
         for j in range(exp_data.shape[1]):
            if exp_data[i,j] < threshold:
               if not (i%mod==mod/2 and j%mod==mod/2):
                  mask[i,j] = 0
   final=total-mask.sum()
   #print('Of {0} pixels, {1} are masked by beamstop and {2} are being skipped for speed.'.format(total,int(starting),int(final-starting)))
   print('Using {0} pixels out of a total of {1}.'.format(total-int(final),total))
   print('{0} are masked by beamstop and {1} are being skipped for speed.'.format(int(starting),int(final-starting)))
   print('Estimated speedup: {0:.3}x.'.format(float(total)/(total-final)))
   return mask

def convert_from_SI():
    '''Copies parameters from SI g.dictionary back to regular g.dictionary.'''
    g.dictionary["radius_1"] = g.dictionary_SI["radius_1"]*10**9
    g.dictionary["radius_2"] = g.dictionary_SI["radius_2"]*10**9
    g.dictionary["z_dim"] = g.dictionary_SI["z_dim"]*10**9
    if g.dictionary["degrees"] == 1: #Converting from radians
       g.dictionary["x_theta"] = g.dictionary_SI["x_theta"]*180/np.pi
       g.dictionary["y_theta"] = g.dictionary_SI["y_theta"]*180/np.pi
       g.dictionary["z_theta"] = g.dictionary_SI["z_theta"]*180/np.pi
    g.dictionary["rho_1"] = g.dictionary_SI["rho_1"]
    g.dictionary["rho_2"] = g.dictionary_SI["rho_2"]

def plot_residuals():
   '''Loads exp data, calculates intensity, and plots the difference [as well as 2 original plots].'''
   plot_all=g.dictionary['plot_fit_tick']
   get_numbers_from_gui()
   load_functions()
   filename = g.dictionary['fit_file']
   if not g.f2py_enabled and not g.opencl_enabled:
      print('Acceleration is NOT enabled!')
   print('{0}: Starting calculation...'.format(time.strftime("%X")))
   exp_data,mask=load_exp_image()
   if g.dictionary['grid_compression'] > 1:
      fast_mask(exp_data,mask,g.dictionary['grid_compression'])
   #calc_intensity=Average_Intensity(mask)
   calc_intensity = normalize(Detector_Intensity(Points_For_Calculation(),mask),mask,True)
   save(calc_intensity,"_calc")     #wrong suffix!!
   err = calc_intensity - exp_data
   ##calc_intensity = normalize(Detector_Intensity(Points_For_Calculation(seed=random_seed),mask),mask,True)
   #err = mask*(exp_data - (calc_intensity + g.dictionary_SI['background']))
   #calc_intensity += g.dictionary_SI'background'] - exp_data  #todo: might be faster
   #calc_intensity *= mask                                     #todo: might be faster
   save(err,"_guess_residuals")
   plot_residuals=np.abs(err)*mask
   print('{0}: Total error = {1:.4}; sum of squares = {2:.4}'.format(time.strftime("%X"),plot_residuals.sum(),np.square(err).sum()))
   if plot_all:
      Fit_plot(exp_data*mask,calc_intensity,plot_residuals)
   else:
      print('Plotting difference.')
      Intensity_plot(plot_residuals,"residuals",'Difference Plot',1)




### End Fitting Functions ###





################          MAKING THE GUI            #################
advanced = ['log_scale2', 'bound2', 'minimum', 'ave_dist', 'scale2',
              'ThreeD2', 'pixels', 'maximum', 'altitude', 'azimuth', 'z_scale']
#All these things will be hiden when advanced/simple options button is pressed

#This refers to 
def hide(event):
   if g.dictionary['advanced'] == 1:
      for x in advanced:
         g.dictionary_in[x].config(state=DISABLED)
      g.dictionary['advanced']=0
      advbutton["text"] = "Advanced Options"
   else:
      for x in advanced:
         g.dictionary_in[x].config(state=NORMAL)
      g.dictionary['advanced']=1
      advbutton["text"] = "Simple Options"


#These are all hidden when Edit/No sequence is pressed.
seq_options = ['s_step', 's_stop', 's_start', 's_var', 'gauss0', 'gauss1', 'SD']
def hide_sequence(event):
   if g.dictionary['seq_hide'] == 1:
      for x in seq_options:
         g.dictionary_in[x].config(state=DISABLED)
      g.dictionary['seq_hide']=0
      seq_button["text"] = "Edit Sequence"
      int_button['text'] = 'Calculate Intensity'
   else:
      for x in seq_options:
         g.dictionary_in[x].config(state=NORMAL)
      g.dictionary['seq_hide']=1
      seq_button["text"] = "No Sequence"
      int_button['text'] = 'Calculate Sequence'

#If you want a number box and a label, use this
def enter_num(variable_name, label, ROW, COL):
    Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
    g.dictionary_in[variable_name] = StringVar()
    g.dictionary_in[variable_name].set(g.dictionary[variable_name])
    g.dictionary_in[variable_name] = Entry(master, textvariable = g.dictionary_in[variable_name])
    if variable_name in advanced and g.dictionary['advanced'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and g.dictionary['seq_hide'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    g.dictionary_in[variable_name].grid(row= ROW, column = COL+1)

#If you want a number box below the label
def enter_vert_num(variable_name, label, ROW, COL):
    Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
    g.dictionary_in[variable_name] = StringVar()
    g.dictionary_in[variable_name].set(g.dictionary[variable_name])
    g.dictionary_in[variable_name] = Entry(master, textvariable = g.dictionary_in[variable_name])
    if variable_name in advanced and g.dictionary['advanced'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and g.dictionary['seq_hide'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    g.dictionary_in[variable_name].grid(row= ROW+1, column = COL)

#If you want a tick box
def tick(variable_name, label, ROW, COL):
    g.dictionary_in[variable_name] = IntVar()
    g.dictionary_in[variable_name].set(int(g.dictionary[variable_name]))
    g.dictionary_in[variable_name+'2'] = Checkbutton(master, text=label, variable=g.dictionary_in[variable_name])
    if variable_name+'2' in advanced and g.dictionary['advanced'] == 0:
       g.dictionary_in[variable_name+'2'].config(state = DISABLED)

    if variable_name+'2' in seq_options and g.dictionary['seq_hide'] == 0:
       g.dictionary_in[variable_name+'2'].config(state = DISABLED)
    g.dictionary_in[variable_name+'2'].grid(row=ROW, column = COL, sticky=W)

#If you want a string entered
def enter_str(variable_name, label, ROW, COL):
    Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
    g.dictionary_in[variable_name] = Entry(master)
    g.dictionary_in[variable_name].insert(0, g.dictionary[variable_name])
    if variable_name in advanced and g.dictionary['advanced'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and g.dictionary['seq_hide'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    g.dictionary_in[variable_name].grid(row= ROW, column = COL + 1)

def enter_text(variable_name, label, WIDTH, HEIGHT, ROW, COL):#For a large textbox
   Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
   g.dictionary_in[variable_name] = Text(master, height = HEIGHT, width = WIDTH)

   g.dictionary_in[variable_name].insert(1.0, g.dictionary[variable_name])
   g.dictionary_in[variable_name].grid(row=ROW+1, column = COL, rowspan = 2, sticky = W)

def radio(variable_name, MODES, ROW, COL): #Radiobutton
   g.dictionary_in[variable_name] = StringVar()
   g.dictionary_in[variable_name].set(g.dictionary[variable_name])
   for name, mode in MODES:
      g.dictionary_in[variable_name+mode] = Radiobutton(master, text=name, variable=g.dictionary_in[variable_name], value=mode)
      g.dictionary_in[variable_name+mode].grid(row=ROW, column = COL, sticky = W)
      if variable_name+mode in advanced and g.dictionary['advanced'] == 0:
          g.dictionary_in[variable_name+mode].config(state = DISABLED)
      if variable_name+mode in seq_options and g.dictionary['seq_hide'] == 0:
          g.dictionary_in[variable_name+mode].config(state = DISABLED)
      COL+=1




###########THE GUI STARTS HERE#####################

if __name__ == "__main__":
   master = Tk()
   master.title("Monte Carlo Small Angle Scattering ({0}), By Max Proft".format(version))

   ### Model Type ###

   ROW = 0
   COL = 0
   Label(master, text = "Model Type", font = "Times 16 bold").grid(row = ROW, column = COL, sticky = W)
   ROW+=1
   
   Label(master, text = "Choose a Monte Carlo Model").grid(row = ROW, column = COL, sticky = W)
   g.dictionary_in['shape'] = StringVar(master)
   g.dictionary_in['shape'].set(MC_num_and_name[g.dictionary['shape']][0])
   OptionMenu(master, g.dictionary_in['shape'], *MC_num_and_name[:,0]).grid(row = ROW, column = COL+1)
   
   ROW+=1
   tick("symmetric", "Radial Symmetry", ROW,COL)
   ROW+=1
   #COL+=1
   tick('Qz',"Small Angle Approx. (Qz=0)", ROW, COL)
   #COL-=1
   
   ROW+=1
   Label(master, text = "Choose an Analytic Model").grid(row = ROW, column = COL, sticky = W)
   g.dictionary_in['analytic'] = StringVar(master)
   g.dictionary_in['analytic'].set(Analytic_options[g.dictionary['analytic']][0])
   OptionMenu(master, g.dictionary_in['analytic'], *Analytic_options[:,0]).grid(row = ROW, column = COL+1)
   
   ### Parameters ###
   
   ROW+=1
   Label(master, text = "Parameters", font = "Times 16 bold").grid(row = ROW, column = COL, sticky = W)
   ROW+=1
   enter_num('num_plots', "Number of Plots to Average", ROW, COL)
   ROW+=1
   enter_num('radius_1', "Radius 1 (nm)", ROW, COL)
   ROW+=1
   enter_num('radius_2', "Radius 2 (nm)", ROW, COL)
   ROW+=1
   enter_num('z_dim', "Length (nm)", ROW, COL)
   ROW+=1
   enter_num('rho_1', 'Rho 1', ROW, COL)
   ROW+=1
   enter_num('rho_2', 'Rho 2', ROW, COL)
   ROW+=1
   MODES_energy_wave = [('Energy (keV)', '0'),('Wavelength (nm)','1'),]
   radio('energy_wavelength_box', MODES_energy_wave, ROW, COL)
   ROW+=1
   enter_num('energy_wavelength', "Energy/Wavelength", ROW, COL)
   ROW+=1
   enter_num('QSize', "Q Range (nm^-1)\nof ENTIRE Detector", ROW, COL)
   ROW+=1
   
   Label(master, text = "Angle of rotation", font = "Times 11 underline").grid(row = ROW, column = COL, sticky = W)
   ROW+=1
   MODES_angle = [('Degrees', '1'),('Radians','0'),]
   radio('degrees', MODES_angle, ROW, COL)
   ROW+=1
   enter_num('x_theta', 'x rotation', ROW, COL)
   ROW+=1
   enter_num('y_theta', 'y rotation', ROW, COL)
   ROW+=1
   enter_num('z_theta', 'z rotation', ROW, COL)
   
   ROW+=3
   enter_num('circ_delta', 'Ring Thickness (pixels)', ROW, COL)
   ROW+=1
   enter_num('theta_delta','Number of Points',ROW, COL)
   ROW+=1
   enter_num('proportional_radius', 'Proportional Radius', ROW, COL)
   ROW+=1
  
   ### Output Options ###
   
   COL+=2
   ROW=0
   Label(master, text="Output Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   COL+=1
   if g.dictionary['advanced'] == 0:
      advbutton = Button(master, text="Advanced Options", font = "Times 12 bold")
   else:
      advbutton = Button(master, text="Simple Options", font = "Times 12 bold")
   advbutton.bind("<Button-1>", hide)
   advbutton.grid(row=ROW, column = COL, pady=4)
   COL-=1
   
   ROW+=1
   tick('save_img', 'Save Images?', ROW, COL)
   COL+=1
   tick('scale', "Make Axes to Scale when\nPlotting Real Space", ROW, COL)
   COL-=1
   ROW+=1
   tick('log_scale', "Plot on a Log Scale?", ROW, COL)
   COL+=1
   tick('ThreeD', "Plot the Intensity in 3D", ROW, COL)
   COL-=1
   ROW+=1
   Label(master, text="Choose Alt/Az for 3D plots (degrees)").grid(row= ROW, column=COL, sticky = W)
   ROW+=1
   enter_num('altitude', "Altitude", ROW, COL)
   ROW+=1
   enter_num('azimuth', "Azimuth", ROW, COL)
   ROW+=1
   enter_num('pixels', "Number of Pixels (x y)", ROW, COL)
   ROW+=1
   enter_num('ave_dist', "Neighbouring Point Distance (nm)", ROW, COL)
   ROW+=1
   enter_num('z_scale','z-direction scaling of\nneighbouring point distance', ROW, COL)
   ROW+=1
   tick('bound', "Upper and Lower Bounds?", ROW, COL)
   g.dictionary_in['bound2']['font'] = "Times 11 underline"
   ROW+=1
   enter_num('minimum', "Minimum (enter it in the form: 3e-7)", ROW, COL)
   ROW+=1
   enter_num('maximum', "Maximum", ROW, COL)
   ROW+=1
   
   ### File Information ###
   
   Label(master, text="File Infomation", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   ROW+=1
   enter_str('title', 'Plot Title', ROW, COL)
   ROW+=1
   enter_str('save_name', 'File Name', ROW, COL)
   ROW+=1
   enter_str('subfolder', 'Subfolder', ROW, COL)
   ROW+=1
   Label(master, text="(No spaces at the start or end of File Name or Subfolder!)").grid(row= ROW, column=COL, columnspan = 2, sticky = W)
   ROW+=1
   ROW+=1
   HEIGHT = 3
   WIDTH = 25
   enter_text("comments", "Description (optional):", WIDTH, HEIGHT, ROW, COL)
   ROW+=3
   
   
   ROW+=1
   
   
   button_row = ROW
   if g.dictionary['seq_hide'] == 0:
      int_button = Button(master, text="Calculate Intensity", command = int_seq, font = "Times 16 bold")
   else:
      int_button = Button(master, text="Calculate Sequence", command = int_seq, font = "Times 16 bold")
   int_button.grid(row=ROW, column = COL,rowspan = 2, pady=4)
   
   COL+=1
   Button(master, text='Replot Intensity', command=view_intensity, font = "Times 16 bold").grid(row=ROW, column = COL,rowspan = 2, pady=4)
   COL-=1
   ROW+=2
   Button(master, text='Real Space', command=plot_points, font = "Times 16 bold").grid(row=ROW, column = COL,rowspan = 2, pady=4)
   COL+=1
   Button(master, text='Plot a Ring', command=circ, font = "Times 16 bold").grid(row=ROW, column = COL, rowspan = 2, pady=4)
   COL-=1
   ROW+=1






   ### Sequence Options ###
   
   COL += 2
   ROW = 0
   Label(master, text="Sequence Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   COL+=1
   if g.dictionary['seq_hide'] == 0:
      seq_button = Button(master, text="Edit Sequence", font = "Times 12 bold")
   else:
      seq_button = Button(master, text="No Sequence", font = "Times 12 bold")
   seq_button.bind("<Button-1>", hide_sequence)
   seq_button.grid(row=ROW, column = COL, pady=4)
   COL-=1
   ROW+=1
   enter_num('s_step', 'Number of Frames', ROW, COL)
   ROW+=1
   enter_str('s_var', "Which Variable?", ROW, COL)
   ROW+=1
   Button(master, text='Common Variables', command=show_sequence_variables).grid(row=ROW, column = COL, sticky=W, pady=4)
   
   ROW+=1
   MODES_gauss = [('Linear Sequence', '0'),('Gaussian', '1'),]
   radio('gauss',MODES_gauss, ROW, COL)
   ROW+=1
   enter_num('s_start', "Sequence Start", ROW, COL)
   ROW+=1
   enter_num('s_stop', "Sequence Stop", ROW, COL)
   ROW+=1
   enter_num('SD', 'Standard Deviation', ROW, COL)
   ROW+=1



   ### Fitting Options ###

   Label(master, text="Fitting Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   ROW += 1
   enter_num('fit_file', "Experimental Data Filename", ROW, COL)
   ROW += 1
   enter_num('center', "Center of Beamstop (x y)", ROW, COL)
   ROW += 1
   enter_num('border', "Additional Cropping", ROW, COL)
   ROW += 1
   enter_num('mask_threshold', "Mask Threshhold", ROW, COL)
   ROW += 1
   COL += 1
   Button(master, text="Plot Exp Data", command = plot_exp_data, font = "Times 16 bold").grid(row=ROW, column=COL, pady=4)
   COL -= 1
   ROW += 1
   Label(master, text="Fit Parameters:").grid(row= ROW, column=COL, columnspan =2, sticky = W)
   ROW += 1
   enter_num('background', "Background Noise", ROW, COL)

   ROW+=1      # These are checkboxes which, if unchecked, will hold fixed fit parameters.
   tick("fit_radius_1", "Radius 1", ROW,COL)
   COL+=1
   tick("fit_radius_2", "Radius 2", ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_rho_1", "Rho 1", ROW,COL)
   COL+=1
   tick("fit_rho_2", "Rho 2", ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_z_dim", "Length", ROW,COL)
   COL+=1
   tick("fit_x_theta", "x rotation", ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_y_theta", "y rotation", ROW,COL)
   COL+=1
   tick("fit_z_theta", "z rotation", ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_background", "background", ROW,COL)
   COL+=1
   tick("fit_other", "unused", ROW,COL)
   COL-=1
   ROW += 1
   ROW += 1
   enter_num('max_iter', "Maximum Iterations (0=default)", ROW, COL)
   ROW += 1
   if g.debug:
      enter_num('update_freq', "Update Interval", ROW, COL)   #TODO: debug
      ROW += 1
   enter_num('grid_compression', "Grid Compression (2, 5, or 10)", ROW, COL)
   ROW += 1
   tick('plot_fit_tick',"Plot Fit Results", ROW, COL)
   COL += 1
   tick('plot_residuals_tick',"Plot Fit Residuals", ROW, COL)
   COL-=1
   ROW += 1
   Button(master, text="Plot Residuals", command = plot_residuals, font = "Times 16 bold").grid(row=ROW, column=COL, pady=4)
   COL+= 1
   Button(master, text="Fit Exp Data", command = perform_fit, font = "Times 16 bold").grid(row=ROW, column=COL, pady=4)
   COL -= 1

   
   mainloop()
   



