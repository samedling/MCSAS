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
      if os.path.isfile('fastmath.so'):
         import fastmath
      elif sys.platform == 'darwin':
         os.system('cp fastmath-OSX10.10_C2DP8700.so fastmath.so')
         import fastmath
      elif sys.platform == 'linux2':
         os.system('cp fastmath-Ubuntu14.10_i7M640.so fastmath.so')
         import fastmath
      print("Accelerating using f2py.")
   except ImportError:
      g.f2py_enabled = False
   #try:
   #   fastmath.module.new_function(<vars>)     #Update this line to check for updates for f2py binary.
   #except AttributeError:
   #   print('Existing f2py binary was out of date; newer version copied.')
   #   vprint('If you ran `make` yourself, run it again for optimal performance.')
   #   if sys.platform == 'darwin':
   #      os.system('cp fastmath-OSX10.10_C2DP8700.so fastmath.so')
   #      import fastmath
   #   elif sys.platform == 'linux2':
   #      os.system('cp fastmath-Ubuntu14.10_i7M640.so fastmath.so')
   #      import fastmath
   #   print("Accelerating using f2py.")
if not g.f2py_enabled and not g.opencl_enabled:
   print("Could not accelerate using either OpenCL or f2py.")
   print("See README for how to install either OpenCL or f2py.")
   print("In the meantime, fitting is not recommended.")

if g.debug:
   np.random.seed([2015])     #Locks random seed to allow for speedtesting.


#These are the default settings
g.dictionary = {'advanced':1, 'altitude':45, 'analytic': 2, 'ave_dist': 0.6, 'azimuth':45, 'bound': 1, 'circ_delta':5, 'comments':'',
              'degrees': 1, 'energy_wavelength': 11, 'energy_wavelength_box': 0, 'gauss':0, 'log_scale': 1, 'maximum': 0.01, 'minimum': 1e-8, 'd_lambda': 2e-4,
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
else:
   root_folder = os.getcwd()  #Mac/Linux
#todo: Check for write access?

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

#from Monte_Carlo_Functions import *
#from Plotting_Functions import *
#from density_formula import *
#from analytic_formula import *

execfile(root_folder+"/Monte_Carlo_Functions.py",globals())
execfile(root_folder+"/Plotting_Functions.py", globals())
execfile(root_folder+"/density_formula.py", globals())
execfile(root_folder+"/analytic_formula.py", globals())
execfile(root_folder+"/fit.py", globals())

#This is the list of all the Monte Carlo Models that you choose from.
##TODO: move this to the same file as the Monte Carlo Models (density_formula.py)
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
                        ["Eliptical Cylinder",15],
                        ["Asym Hex Pyramid",16]
                        ])
MC_num_and_name_dict = {x[0]:x[1] for x in MC_num_and_name} #This is needed, so that when an option is chosen, we can find the shape number.

#This is the list of all the analytic models that you choose from.
##TODO: move this to the same file as the Analytic Models (analytic_formula.py)
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
          try:
             g.dictionary[x] = g.dictionary_in[x].get(1.0,END).rstrip()
          except AttributeError: #'int' object has no attribute 'get'...for when testing and there's no actual GUI
             g.dictionary[x] = g.dictionary_in[x]
       elif x=='advanced' or x== 'seq_hide':
          g.dictionary[x] = g.dictionary[x]
       elif x=='shape':
          try:
             g.dictionary[x] = MC_num_and_name_dict[g.dictionary_in['shape'].get()]
          except AttributeError: #'int' object has no attribute 'get'...for when testing and there's no actual GUI
             g.dictionary[x] = g.dictionary_in['shape']
       elif x=='analytic':
          try:
             g.dictionary[x] = Analytic_dict[g.dictionary_in['analytic'].get()]
          except AttributeError: #'int' object has no attribute 'get'...for when testing and there's no actual GUI
             g.dictionary[x] = g.dictionary_in['analytic']
       else:
          try:
             g.dictionary[x] = g.dictionary_in[x].get() 
          except AttributeError:
             g.dprint("{0} could not be imported from GUI.".format(x))
    
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
   #ROW+=1
   #tick('Qz',"Small Angle Approx. (Qz=0)", ROW, COL)
   COL+=1
   tick('Qz',"Small Angle Approx.", ROW, COL)
   COL-=1
   
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
   enter_num('d_lambda', "Wavelength Spread (2e-4)", ROW, COL)
   ROW+=1
   enter_num('QSize', "Detector Q Range (nm^-1)", ROW, COL)
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
   



