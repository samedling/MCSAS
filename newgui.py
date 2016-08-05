#!/usr/bin/python
version = '0.6.0'
updated = '5 Aug 2016'

print('Starting MCSAS v{0} (updated {1}).'.format(version,updated))


try:
   from Tkinter import *
except ImportError:
   from tkinter import *

import os, sys, pickle, random, pylab, time
import numpy as np
from scipy.special import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

try:
   from PIL import Image
except ImportError:
   import Image   #not sure if this will work
   #from scipy import misc     #alternative to PIL?

from scipy.optimize import leastsq
#from scipy import ndimage  #possible smoothing of exp_data before viewing

import global_vars as g

if sys.version_info[0] > 2:
   print("This program has not been tested with Python 3. Python 2.7 is recommended.")
elif sys.version_info[1] < 7:
   print("This program requires Python version 2.7 or later to run. Exiting.")
   sys.exit()

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
         print('No fortran binary found; copying mac binary.')
         import fastmath
      elif sys.platform == 'linux2':
         os.system('cp fastmath-Ubuntu14.04_i7-4770.so fastmath.so')
         #os.system('cp fastmath-Ubuntu14.10_i7M640.so fastmath.so')
         print('No fortran binary found; copying linux binary.')
         import fastmath
      print("Accelerating using f2py.")
   except ImportError:
      g.f2py_enabled = False
   try:
     if fastmath.version.number() < 0.405:     #Update this to check f2py binary updates. Value should be 0.005 less than fortran.
        print('Existing f2py binary was out of date; newer version copied.')
        g.vprint('If you ran `make` yourself, run it again for optimal performance.')
        if sys.platform == 'darwin':
           os.system('cp fastmath-OSX10.10_C2DP8700.so fastmath.so')
           import fastmath
        elif sys.platform == 'linux2':
           os.system('cp fastmath-Ubuntu14.10_i7M640.so fastmath.so')
           import fastmath
        print("Accelerating using f2py.")
   except AttributeError:
     print('Existing f2py binary was out of date; newer version copied.')
     g.vprint('If you ran `make` yourself, run it again for optimal performance.')
     if sys.platform == 'darwin':
        os.system('cp fastmath-OSX10.10_C2DP8700.so fastmath.so')
        import fastmath
     elif sys.platform == 'linux2':
        os.system('cp fastmath-Ubuntu14.10_i7M640.so fastmath.so')
        import fastmath
     print("Accelerating using f2py.")
   except NameError: #For when fortran binary isn't compatible.
      print('Fortran binary is not compatible. Please run make to compile.')
      g.f2py_enabled = False
if not g.f2py_enabled and not g.opencl_enabled:
   print("Could not accelerate using either OpenCL or f2py.")
   print("See README for how to install either OpenCL or f2py.")
   print("In the meantime, fitting is not recommended.")

if g.debug:
   np.random.seed([2015])     #Locks random seed to allow for speedtesting.


#These are the default settings
g.dictionary = {'advanced':1, 'altitude':45, 'analytic': 2, 'ave_dist': 1.0, 'azimuth':45, 'bound': 1, 'circ_delta':5, 'comments':'',
              'degrees': 1, 'energy_wavelength': 11, 'energy_wavelength_box': 0, 'gauss':0, 'log_scale': 1, 'maximum': 0.01, 'minimum': 1e-8, 'd_lambda': 2e-4,
              'num_plots': 1, 'pixels': (200,200), 'proportional_radius':0.5, 'QSize': 6,'Qz': 0, 'radius_1': 5.0, 'radius_2': 2.5, 'rho_1': 1.0, 'rho_2': -0.5,
              'save_img':1, 'save_name': 'save_name', 'scale': 1,'SD':1, 'seq_hide':1, 'shape': 2, 's_start': 0, 's_step': 2,
              's_stop': 1, 'subfolder':'subfolder', 's_var': 'x_theta', 'symmetric': 0, 'num':1, 'length_2':0,
              'theta_delta':20, 'ThreeD': 0, 'title': 'title', 'x_theta': 0,'y_theta': 0,'z_theta': 0,'z_dim': 100,'z_scale':1,#}
              'fit_file': 'fit_file', 'center': (0,0), 'border': 0, 'max_iter': 1000, 'update_freq': 0, 'plot_fit_tick': 1, 'plot_residuals_tick': 1, 'mask_threshold': 10, 'background': 2e-5, 'grid_compression': 0,
              'fit_radius_1': 1, 'fit_radius_2': 0, 'fit_rho_1': 1, 'fit_rho_2': 0, 'fit_z_dim': 1, 'fit_x_theta': 1, 'fit_y_theta': 1, 'fit_z_theta': 1, 'fit_background': 1, 'fit_num': 0, 'fit_length_2':0,
              'xinter':100,'yinter':100,'numinter':10,'save_points':0,
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
    for x in g.dictionary:
      if x not in d: 
         d[x]=g.dictionary[x]
    g.dictionary = d
except:
    print("Previously used variables could not be loaded. \nUsing default settings instead.")
    with open(root_folder+"/default.txt", 'wb') as f:
       pickle.dump(g.dictionary, f)

g.dictionary = {x:g.dictionary[x] for x in g.dictionary} #This contains the unaltered parameters
g.dictionary_in = {a:g.dictionary[x] for x in g.dictionary for a in [x, x+'2']} #This contains raw data from the GUI. The +'2' is so that i can control checkboxes
g.dictionary_SI = {x:g.dictionary[x] for x in g.dictionary} #this g.dictionary has the parameters after they have been converted to SI units

#from Monte_Carlo_Functions import *
#from Plotting_Functions import *
#from density_formula import *
#from analytic_formula import *

try:
   execfile(root_folder+"/Monte_Carlo_Functions.py",globals())
   execfile(root_folder+"/Plotting_Functions.py", globals())
   execfile(root_folder+"/density_formula.py", globals())
   execfile(root_folder+"/analytic_formula.py", globals())
   execfile(root_folder+"/fit.py", globals())
except:
   exec(open(root_folder+"/Monte_Carlo_Functions.py").read(),globals())
   exec(open(root_folder+"/Plotting_Functions.py").read(),globals())
   exec(open(root_folder+"/density_formula.py").read(),globals())
   exec(open(root_folder+"/analytic_formula.py").read(),globals())
   exec(open(root_folder+"/fit.py").read(),globals())

#This is the list of all the analytic models that you choose from.
##TODO: move this to the same file as the Analytic Models (analytic_formula.py)
Analytic_options = np.array([["None",0],
                        ["Sphere",1],
                        ["Cylinder",2],
                        ["Core shell cylinder",3],
                        ["Gaussian",4]
                        ])
Analytic_dict = {x[0]:x[1] for x in Analytic_options} #This is needed, so that when an option is chosen, we can find the shape number.


def xy_dim(): #Since x_dim and y_dim are not defined by the user any more, the function to define them is here.
   if g.model_parameters[g.dictionary['shape']][2][1] == 'unused':
      g.dictionary_SI['x_dim'] = 2.*g.dictionary_SI['radius_1']
      g.dictionary_SI['y_dim'] = 2.*g.dictionary_SI['radius_1']
   else:
      g.dictionary_SI['x_dim'] = 2.*max(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'])
      g.dictionary_SI['y_dim'] = 2.*max(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'])
   
   

def get_numbers_from_gui():
    #Here I get all parameters from the GUI and put them into g.dictionary
    for x in g.dictionary:
       if x=='comments':
          try:
             g.dictionary[x] = g.dictionary_in[x].get(1.0,END).rstrip()
          except AttributeError: #'int' object has no attribute 'get'...for when testing and there's no actual GUI
             #g.dictionary[x] = g.dictionary_in[x]
             pass
       elif x=='advanced' or x== 'seq_hide':
          g.dictionary[x] = g.dictionary[x]
       elif x=='shape':
          try:
             g.dictionary[x] = g.MC_num_and_name_dict[g.dictionary_in['shape'].get()]
          except AttributeError: #'int' object has no attribute 'get'...for when testing and there's no actual GUI
             #g.dictionary[x] = g.dictionary_in['shape']
             pass
       elif x=='analytic':
          try:
             g.dictionary[x] = Analytic_dict[g.dictionary_in['analytic'].get()]
          except AttributeError: #'int' object has no attribute 'get'...for when testing and there's no actual GUI
             #g.dictionary[x] = g.dictionary_in['analytic']
             pass
       else:
          try:
             g.dictionary[x] = g.dictionary_in[x].get() 
          except (TclError, AttributeError):#If it can't load it, if you open then close a window for example, it just uses the existing value.
             #g.dprint("{0} could not be imported from GUI.".format(x))
             pass
    
    for x in g.dictionary:
        try:
            g.dictionary[x] = float(g.dictionary[x]) #I am turning the numbers from strings to floats.
            if g.dictionary[x]==int(g.dictionary[x]):#I am turning floats to integers if they are the same
               g.dictionary[x] = int(g.dictionary[x])
        except:
            None
    if not os.path.exists(root_folder+'/'+g.dictionary['subfolder']):#making the subfolder, if it doesn't exist
       os.makedirs(root_folder+'/'+g.dictionary['subfolder'])
       #time.sleep(2) #Making the subfolder takes a few seconds, so we need to delay the program, otherwise it will try save things into the folder before it is made.

    make_SI_dict()


def make_SI_dict():
    g.dictionary_SI = {x: g.dictionary[x] for x in g.dictionary}
    g.dictionary_SI['path_to_subfolder'] = os.path.join(root_folder,g.dictionary['subfolder'],g.dictionary['save_name']) #This is for convienience

    #Converting to SI units.
    g.dictionary_SI["z_dim"] = g.dictionary["z_dim"]*10**-9
    g.dictionary_SI["length_2"] = g.dictionary["length_2"]*10**-9
    g.dictionary_SI["ave_dist"] = g.dictionary["ave_dist"]*10**-9
    g.dictionary_SI["travel"] = g.dictionary_SI["ave_dist"]
    g.dictionary_SI["radius_1"] = g.dictionary["radius_1"]*10**-9
    g.dictionary_SI["radius_2"] = g.dictionary["radius_2"]*10**-9
    xy_dim()#defining x_dim and y_dim - dependent of radius_1
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
   print( g.dictionary_SI['path_to_subfolder'])
   print( g.dictionary['comments'])
   print( extra)
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
                ["num","Number of Pieces"],
                ["length_2","Other Length"],
                ["pixels","pixels on the detector"],
                ["QSize","Q Range of the entire detector"],
                ["radius_1","Radius 1"],
                ["radius_2","Radius 2"],
                ["rho_1","Rho 1"],
                ["rho_2","Rho 2"],
                ["x_theta","Rotation in x direction"],
                ["y_theta","Rotation in y direction"],
                ["z_dim","Length"],
                ['z_scale', 'Z-direction Scaling of\nNeighbouring Point Distance'],
                ["z_theta","Rotation in z Direction"],
                ]
ALLVAR_LIST= []
ALLVAR_NAMES= []
for i in range(len(ALLVARIABLES)):
   ALLVAR_LIST.append(ALLVARIABLES[i][0])
   ALLVAR_NAMES.append(ALLVARIABLES[i][1])

def show_sequence_variables(): #Common Variables button, displays ALLVARIABLES, above.
   get_numbers_from_gui()
   #rename_parameters()
   for i in range(len(g.var_list)):
      ALLVAR_NAMES[ALLVAR_LIST.index(g.var_list[i])] = g.var_names[i]
   dens_options = Tk()
   dens_options.title("Common Variables")
   #Label(dens_options, text = "Choose a Variable to Change", font = "Times 14 bold underline").grid(row = 0, column = 0, columnspan = 2, sticky = W)
   Label(dens_options, text = "Variable", font = "Times 11 underline").grid(row = 0, column = 0, sticky = W)
   Label(dens_options, text = "Description", font = "Times 11 underline").grid(row = 0, column = 1, sticky = W)
   for ROW in range(len(ALLVARIABLES)):
      Label(dens_options, text = ALLVAR_LIST[ROW]).grid(row = ROW+2, column = 0, sticky = W)
      Label(dens_options, text = ALLVAR_NAMES[ROW]).grid(row = ROW+2, column = 1, sticky = W)
   dens_options.mainloop()

def plot_points(): #This runs the Real Space to plot the points in Real Space
    get_numbers_from_gui()
    save_vars_to_file("Plot Points")
    if g.dictionary['seq_hide'] == 1:
       if g.dictionary['gauss']==0:#if 'gauss', it takes a random number from the gaussian distribution.
          current_value = g.dictionary['s_start']
       else:
          current_value = np.random.normal(loc = g.dictionary[g.dictionary['s_var']], scale = g.dictionary['SD'])
       change_units(current_value)
    
    Points = Points_For_Calculation()
    if g.dictionary_SI['save_points']==1:
          g.vprint("Saving Points")
          save(Points,"Points")
    Points_Plot(Points, 'points', 1)
    clear_mem()
    print("Program Finished")

def view_intensity(): #This allows you to view a premade intensity
    get_numbers_from_gui()
    radial_intensity = pylab.loadtxt(g.dictionary_SI['path_to_subfolder']+"radial_intensity.csv", delimiter=",")
    radial_intensity_plot(radial_intensity, "radial", g.dictionary_SI['title'], 0)
    Intensity = pylab.loadtxt(g.dictionary_SI['path_to_subfolder']+"intensity.csv", delimiter=",")
    Intensity_plot(Intensity, "intensity", g.dictionary_SI['title'], 1)
    clear_mem()
    print( "Program Finished")

def calc_coherence_length():
   get_numbers_from_gui()
   if g.dictionary['d_lambda']:
      coherence_length = 2*np.pi/(g.dictionary_SI['EHC']*g.dictionary['d_lambda'])
      print('Coherence length is {0:6.4} nm.'.format(coherence_length*10**9))
   else:
      coherence_length = 1
      print('Coherence length is {0} m.'.format(coherence_length))


def slow_intensity(): #This makes an intensity the accuate, slow way.
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file("Monte Carlo Intensity")
    print( "START TIME: "+time.strftime("%X"))
    Intensity = normalize(Accurate_Intensity(Points_For_Calculation(sort=1)))
    print( "END TIME: "+time.strftime("%X"))
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    if g.dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    
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
    print( "Program Finished")
    
def interparticle():
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file('Interparticle Scattering')
    Intensity = inter_intensity()
    print( "END TIME: "+time.strftime("%X"))
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    if g.dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
 

def sequence(): #This makes a sequence of intensities
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file("Monte Carlo Sequence")
    for frame_num in range(int(g.dictionary['s_step'])):
       sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
       sim_info.write("\nFrame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step'])))
       sim_info.close()

       print( "\nmaking frame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step'])))
       if g.dictionary['gauss']==0:
          try:
             current_value = (g.dictionary['s_stop']-g.dictionary['s_start'])*frame_num/(1.*g.dictionary['s_step']-1.)+1.*g.dictionary['s_start']
          except ZeroDivisionError:
             current_value = g.dictionary['s_start']
       else:
          current_value = np.random.normal(loc = g.dictionary[g.dictionary['s_var']], scale = g.dictionary['SD'])
       change_units(current_value)
       Intensity = Average_Intensity(seqnum=str(frame_num+1))
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
    print( "Program Finished")


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
   print( "Program Finished")



def theory_seq(): #This plots a sequence created with the analytic model
    global sim_info
    get_numbers_from_gui()
    save_vars_to_file("Analytic Sequence")
    for frame_num in range(int(g.dictionary['s_step'])):
       sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
       sim_info.write("\nFrame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step'])))
       sim_info.close()

       print( "\nmaking frame " + str(frame_num+1) + " of " + str(int(g.dictionary['s_step'])))
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
    print( "Program Finished")



def circ(): #This plots a the angle at a fixed radius
   get_numbers_from_gui()
   Intensity = np.asarray(pylab.loadtxt(g.dictionary_SI['path_to_subfolder']+"intensity.csv", delimiter=","))
   data = plotting_circle(Intensity)
   radial_intensity_plot(data, "theta"+str(g.dictionary['radius_2']), g.dictionary['title']+" "+str(g.dictionary['radius_2']), 0)
   angle_plot(data, "Angle"+str(g.dictionary['radius_2']), g.dictionary['title']+" "+str(g.dictionary['radius_2']), 1)
   print( "finshied")
   
def calc_int():
   '''Calculates intensity (average if # plots > 1).'''
   get_numbers_from_gui() #otherwise it uses the previous 'shape' number
   if g.dictionary['shape'] ==0:
      theory_plot()
   else:
      make_intensity()

def calc_seq():
   '''Calculates intensity for sequence (averaging each time if # plots > 1.'''
   get_numbers_from_gui() #otherwise it uses the previous 'shape' number
   if g.dictionary['s_var'] in g.dictionary:
      if g.dictionary['shape']==0:
         theory_seq()
      else:
         sequence()
   else:
      print('No such variable {0}. Check var list.'.format(g.dictionary['s_var']))

#def int_seq(): #This is the button, it runs a sequence or a single image depending on whether or not you can edit a sequence (For both analytic models and Monte Carlo Models)
   #if g.dictionary['seq_hide'] == 0:
      #if g.dictionary['shape'] ==0:
         #theory_plot()
      #else:
         #make_intensity()
   #else:
      #if g.dictionary['shape']==0:
         #theory_seq()
      #else:
         #sequence()





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
def enter_num(variable_name, label, window, ROW, COL):
    g.labels[variable_name] = Label(window, text=label)
    g.labels[variable_name].grid(row= ROW, column=COL, sticky = W)
    temp = StringVar()
    g.dictionary_in[variable_name] = Entry(window, textvariable = temp)
    temp.set(g.dictionary[variable_name])
    #g.dictionary_in[variable_name] = StringVar()
    #g.dictionary_in[variable_name].set(g.dictionary[variable_name])
    #g.dictionary_in[variable_name] = Entry(window, textvariable = g.dictionary_in[variable_name])
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
def tick(variable_name, label, window, ROW, COL,columnspan=1):
    g.dictionary_in[variable_name] = IntVar()
    g.dictionary_in[variable_name].set(int(g.dictionary[variable_name]))
    g.dictionary_in[variable_name+'2'] = Checkbutton(window, text=label, variable=g.dictionary_in[variable_name])
    if variable_name+'2' in advanced and g.dictionary['advanced'] == 0:
       g.dictionary_in[variable_name+'2'].config(state = DISABLED)

    if variable_name+'2' in seq_options and g.dictionary['seq_hide'] == 0:
       g.dictionary_in[variable_name+'2'].config(state = DISABLED)
    g.dictionary_in[variable_name+'2'].grid(row=ROW, column = COL, columnspan=2,sticky=W)

#If you want a string entered
def enter_str(variable_name, label, window, ROW, COL):
    Label(window, text=label).grid(row= ROW, column=COL, sticky = W)
    g.dictionary_in[variable_name] = Entry(window)
    g.dictionary_in[variable_name].insert(0, g.dictionary[variable_name])
    if variable_name in advanced and g.dictionary['advanced'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and g.dictionary['seq_hide'] == 0:
       g.dictionary_in[variable_name].config(state = DISABLED)
    g.dictionary_in[variable_name].grid(row= ROW, column = COL + 1)

def enter_text(variable_name, label, WIDTH, HEIGHT, ROW, COL, columns=2):#For a large textbox
   Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
   g.dictionary_in[variable_name] = Text(master, height = HEIGHT, width = WIDTH)

   g.dictionary_in[variable_name].insert(1.0, g.dictionary[variable_name])
   g.dictionary_in[variable_name].grid(row=ROW+1, column = COL, rowspan = 2, columnspan = columns, sticky = W)

def radio(variable_name, MODES, window, ROW, COL): #Radiobutton
   g.dictionary_in[variable_name] = StringVar()
   g.dictionary_in[variable_name].set(g.dictionary[variable_name])
   for name, mode in MODES:
      g.dictionary_in[variable_name+mode] = Radiobutton(window, text=name, variable=g.dictionary_in[variable_name], value=mode)
      g.dictionary_in[variable_name+mode].grid(row=ROW, column = COL, sticky = W)
      if variable_name+mode in advanced and g.dictionary['advanced'] == 0:
          g.dictionary_in[variable_name+mode].config(state = DISABLED)
      if variable_name+mode in seq_options and g.dictionary['seq_hide'] == 0:
          g.dictionary_in[variable_name+mode].config(state = DISABLED)
      COL+=1



def save_fitparam():
   get_num_from_gui()

def detector_parameters():
    global det_window,master
    det_window = Toplevel(master)
    det_window.title("Detector Options")
    new_frame = Frame(det_window)
    ROW = 0
    COL = 0
    
    MODES_energy_wave = [('Energy (keV)', '0'),('Wavelength (nm)','1'),]
    radio('energy_wavelength_box', MODES_energy_wave, det_window, ROW, COL)
    ROW+=1
    enter_num('energy_wavelength', "Energy/Wavelength", det_window, ROW, COL)
    ROW+=1
    enter_num('QSize', "Detector Q Range (nm^-1)", det_window, ROW, COL)
    ROW+=1
    Label(det_window, text="Coherence Length", font = "Times 16 bold").grid(row= ROW, column=COL, columnspan=2, sticky = W) 
    ROW+=1
    enter_num('d_lambda', "Wavelength Spread (2e-4)", det_window, ROW, COL)
    ROW+=1
    Button(det_window, text='Print Coherence Length', command=calc_coherence_length, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)
    COL+=1
    Button(det_window, text='SLOW Calculate', command=slow_intensity, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)
    COL-=1



def model_parameters():
    global model_window,master
    model_window = Toplevel(master)
    model_window.title("Parameter Options")
    new_frame = Frame(model_window)
    ROW = 0
    COL = 0
    shape = g.dictionary_in['shape'].get()
    
    if shape:
        return
    
    
def rename_parameters(event):
    global master,parameter_start_row
    ROW = parameter_start_row
    COL = 0
    
    g.var_list=['radius_1','radius_2','z_dim','rho_1','rho_2','num','length_2']
    #print g.dictionary_in['shape'].get()
    shape = np.int(g.MC_num_and_name_dict[g.dictionary_in['shape'].get()])
    g.var_names=g.model_parameters[shape][2]
    for i in range(len(g.var_list)):
        g.labels[g.var_list[i]].config(text=g.var_names[i])
        if g.var_names[i] == "unused":
           g.dprint("{0} is unused.".format(g.var_list[i]))
           g.dictionary_in[g.var_list[i]].config(state=DISABLED)
           g.labels[g.var_list[i]].config(state=DISABLED)
           try:
               g.dictionary_in['fit_'+g.var_list[i]+'2'].config(state=DISABLED)
               g.dictionary_in['fit_'+g.var_list[i]+'2'].config(text=g.var_names[i])
           except AttributeError:
               pass
        else:
           g.dictionary_in[g.var_list[i]].config(state=NORMAL)
           g.labels[g.var_list[i]].config(state=NORMAL)
           try:
               g.dictionary_in['fit_'+g.var_list[i]+'2'].config(state=NORMAL)
               g.dictionary_in['fit_'+g.var_list[i]+'2'].config(text=g.var_names[i])
           except AttributeError:
               pass


def sequence_parameters():
   
    global seq_window,master,seq_button
    seq_window = Toplevel(master)
    seq_window.title("Sequence Options")
    new_frame = Frame(seq_window)
    ROW = 0
    COL = 0
   
    # Label(seq_window, text="Sequence Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
    #if g.dictionary['seq_hide'] == 0:
    #   seq_button = Button(seq_window, text="Edit Sequence", font = "Times 12 bold")
    #else:
    #   seq_button = Button(seq_window, text="No Sequence", font = "Times 12 bold")
    #seq_button.bind("<Button-1>", hide_sequence)
    #seq_button.grid(row=ROW, column = COL, pady=2)
    #COL+=1
    #Button(seq_window, text='Common Variables', command=show_sequence_variables).grid(row=ROW, column = COL, sticky=W, pady=2)    
    #COL-=1
    #ROW+=1
    enter_str('s_var', "Which Variable?", seq_window, ROW, COL)
    ROW+=1
    enter_num('s_step', 'Number of Frames', seq_window, ROW, COL)
    ROW+=1
   
    ROW+=1
    MODES_gauss = [('Linear Sequence', '0'),('Gaussian', '1'),]
    radio('gauss',MODES_gauss, seq_window, ROW, COL)
    ROW+=1
    enter_num('s_start', "Sequence Start", seq_window, ROW, COL)
    ROW+=1
    enter_num('s_stop', "Sequence Stop", seq_window, ROW, COL)
    ROW+=1
    enter_num('SD', 'Standard Deviation', seq_window, ROW, COL)
    ROW+=1
    #Button(seq_window, text='Common Variables', command=show_sequence_variables).grid(row=ROW, column = COL, sticky=W, pady=2)    
    #calc_button.grid(row=ROW, column = COL, pady=2, columnspan=2)
    Button(seq_window, text='Common\nVariables', command=show_sequence_variables).grid(row=ROW, column = COL, pady=2)    
    COL+=1
    calc_button = Button(seq_window, text="Calculate\nSequence", command = calc_seq, font = "Times 16 bold")
    calc_button.grid(row=ROW, column = COL, pady=2)
    ROW+=1


def output_options():
   global output_window,master
   output_window = Toplevel(master)
   output_window.title("Output Options")
   new_frame = Frame(output_window)
   ROW=0
   COL=0

   tick('save_points', "Save Points Used in Calculation", output_window,ROW,COL,2)
   ROW+=1
   tick('scale', "Scale Axes When Plotting Real Space", output_window, ROW, COL, 2)
   ROW+=1
   tick('log_scale', "Plot on a Log Scale?", output_window, ROW, COL, 2)
   ROW+=1
   tick('ThreeD', "Plot the Intensity in 3D", output_window, ROW, COL, 2)
   ROW+=1
   Label(output_window, text="Choose Alt/Az for 3D plots (degrees)").grid(row= ROW, column=COL, columnspan=2, sticky = W)
   ROW+=1
   enter_num('altitude', "Altitude", output_window, ROW, COL)
   ROW+=1
   enter_num('azimuth', "Azimuth", output_window, ROW, COL)
   ROW+=1
   #enter_num('ave_dist', "Neighbouring Point Distance (nm)", output_window, ROW, COL)
   #ROW+=1
   #enter_num('z_scale','z-direction scaling of\nneighbouring point distance', output_window, ROW, COL)
   #ROW+=1

    
def ring_options():
    global ring_window,master
    ring_window = Toplevel(master)
    ring_window.title("Ring Options")
    new_frame = Frame(ring_window)
    ROW = 0
    COL = 0

    enter_num('circ_delta', 'Ring Thickness (pixels)', ring_window, ROW, COL)
    ROW+=1
    enter_num('theta_delta','Number of Points',ring_window, ROW, COL)
    ROW+=1
    enter_num('proportional_radius', 'Proportional Radius', ring_window, ROW, COL)
    ROW+=1

    Button(ring_window, text='Plot a Ring', command=circ, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)

def interparticle_options():
    global inter_window, master
    inter_window = Toplevel(master)
    inter_window.title("Inter-Particle Scattering")
    new_frame = Frame(inter_window)
    ROW=0
    COL=0
    
    enter_num('xinter','Box X Dimension (nm)',inter_window,ROW,COL)
    ROW+=1
    enter_num('yinter','Box Y Dimension (nm)',inter_window,ROW,COL)
    ROW+=1
    enter_num('numinter','Number of Particles',inter_window,ROW,COL)
    ROW+=1
    Button(inter_window,text='Calculate', command=interparticle, font = "Times 16 bold").grid(row=ROW,column=COL,sticky=W,pady=2)

def run_file():
   global test_file
   filename=test_file.get()
   if filename:
      execfile(root_folder+"/"+filename,globals())
   else:
      execfile(root_folder+"/test.py",globals())


def run_code():
   global test_code
   print('Printing code:')
   print(test_code.get(1.0,END).rstrip())
   print('Running code...')
   exec( test_code.get(1.0,END).rstrip())
   print('Done.')


def open_debug():
   global debug_window,master
   global test_file,test_code
   debug_window = Toplevel(master)
   debug_window.title("Debug Window")
   new_frame = Frame(debug_window)
   ROW = 0
   COL = 0

   #Button(debug_window, text='Run test.py', command=run_file, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)
   #ROW+=1

   Label(debug_window, text='Filename').grid(row= ROW, column=COL, sticky = W)
   test_file = Entry(debug_window)
   test_file.insert(0, 'test.py')
   test_file.grid(row= ROW, column = COL + 1)
   COL+=2

   Button(debug_window, text='Run File', command=run_file, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)
   ROW+=1
   COL-=2
   
   Label(debug_window, text='Arbitrary Python Code').grid(row= ROW, column=COL, columnspan = 2, sticky = W)
   test_code = Text(debug_window, height = 5, width = 60)
   #test_code.insert(1.0,'#insert code here')
   test_code.insert(1.0,'#get_numbers_from_gui()\n#insert code here\n')
   test_code.grid(row=ROW+1, column = COL, rowspan = 2, columnspan = 3, sticky = W)
   COL+=2

   #enter_text.get(1.0,END).rstrip()
   Button(debug_window, text='Run Python Code', command=run_code, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)
   ROW+=2




###########THE GUI STARTS HERE#####################

if __name__ == "__main__":
   #master = Tk()
   #master.title("Monte Carlo Small Angle Scattering ({0}), By Max Proft".format(version))
   master = Tk()
   master.title("Monte Carlo Small Angle Scattering ({0})".format(version))
   root = Frame(master)

   ### Model Type ###

   ROW = 0
   COL = 0
   Label(master, text = "Model Type", font = "Times 16 bold").grid(row = ROW, column = COL, sticky = W)
   ROW+=1
   tick('Qz',"Small Angle Approx.", master, ROW, COL)
   COL+=1
   tick("symmetric", "Radial Symmetry", master,ROW,COL)
   COL-=1



   ROW+=1
   
   Label(master, text = "Monte Carlo Model: ").grid(row = ROW, column = COL, sticky = W)
   g.dictionary_in['shape'] = StringVar(master)
   g.dictionary_in['shape'].set(g.MC_num_and_name[g.dictionary['shape']][0])
   OptionMenu(master, g.dictionary_in['shape'], *g.MC_num_and_name[:,0]).grid(row = ROW, column = COL+1)
   #MCModel = OptionMenu(master, g.dictionary_in['shape'], *g.MC_num_and_name[:,0])
   #MCModel.grid(row = ROW, column = COL+1)
   #MCModel.bind("<Button-2>",rename_parameters)
   
   
   ROW+=1
   Label(master, text = "Analytic Model: ").grid(row = ROW, column = COL, sticky = W)
   g.dictionary_in['analytic'] = StringVar(master)
   g.dictionary_in['analytic'].set(Analytic_options[g.dictionary['analytic']][0])
   OptionMenu(master, g.dictionary_in['analytic'], *Analytic_options[:,0]).grid(row = ROW, column = COL+1)
   
   ### Parameters ###
   

   ROW+=1
   Label(master, text = "Model Parameters", font = "Times 16 bold").grid(row = ROW, column = COL, sticky = W)
   COL+=1
   parameter_help = Button(master, text='Update Parameters', font = "Times 12 bold")
   parameter_help.bind("<Button-1>", rename_parameters)
   parameter_help.grid(row=ROW, column = COL, pady=2)   
   COL-=1
   ROW+=1
   parameter_start_row=ROW
   enter_num('radius_1', "Radius 1 (nm)", master, ROW, COL)
   ROW+=1
   enter_num('radius_2', "Radius 2 (nm)", master, ROW, COL)
   ROW+=1
   enter_num('z_dim', "Length (nm)", master, ROW, COL)
   ROW+=1
   enter_num('rho_1', 'Rho 1', master, ROW, COL)
   ROW+=1
   enter_num('rho_2', 'Rho 2', master, ROW, COL)
   ROW+=1

   enter_num('num', 'Number', master, ROW, COL)
   ROW+=1
   enter_num('length_2', 'Length 2 (nm)', master, ROW, COL)
   ROW+=1
   
   Label(master, text = "Rotation Parameters", font = "Times 16 bold").grid(row = ROW, column = COL, sticky = W)
   ROW+=1
   MODES_angle = [('Degrees', '1'),('Radians','0'),]
   radio('degrees', MODES_angle, master, ROW, COL)
   ROW+=1
   enter_num('x_theta', 'x rotation', master, ROW, COL)
   ROW+=1
   enter_num('y_theta', 'y rotation', master, ROW, COL)
   ROW+=1
   enter_num('z_theta', 'z rotation', master, ROW, COL)
   
   # ROW+=1
   # MODES_energy_wave = [('Energy (keV)', '0'),('Wavelength (nm)','1'),]
   # radio('energy_wavelength_box', MODES_energy_wave, master, ROW, COL)
   # ROW+=1
   # enter_num('energy_wavelength', "Energy/Wavelength", master, ROW, COL)
   # ROW+=1
   # enter_num('d_lambda', "Wavelength Spread (2e-4)", master, ROW, COL)
   # ROW+=1
   # enter_num('QSize', "Detector Q Range (nm^-1)", master, ROW, COL)
   
   # ROW+=3
   # enter_num('circ_delta', 'Ring Thickness (pixels)', master, ROW, COL)
   # ROW+=1
   # enter_num('theta_delta','Number of Points',master, ROW, COL)
   # ROW+=1
   # enter_num('proportional_radius', 'Proportional Radius', master, ROW, COL)
   # ROW+=1
  
   ### Output Options ###
   
   COL+=2
   ROW=0
   Label(master, text="Output Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   COL+=1
   tick('save_img', 'Save Images?', master, ROW, COL)
   ROW+=1
   COL-=1

   ### File Information ###
   
   #Label(master, text="File Infomation", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   #ROW+=1
   enter_str('title', 'Plot Title', master, ROW, COL)
   ROW+=1
   enter_str('save_name', 'File Name', master, ROW, COL)
   ROW+=1
   enter_str('subfolder', 'Subfolder', master, ROW, COL)
   ROW+=1
   # Label(master, text="(No spaces at the start or end of File Name or Subfolder!)").grid(row= ROW, column=COL, columnspan = 2, sticky = W)
   #ROW+=1
   HEIGHT = 3
   WIDTH = 50
   enter_text("comments", "Description (optional):", WIDTH, HEIGHT, ROW, COL, 2)
   ROW+=3
   
   enter_num('num_plots', "Number of Plots to Average", master, ROW, COL)
   ROW+=1

   #tick('scale', "Scale Axes When\nPlotting Real Space", master, ROW, COL)
   #ROW+=1
   #tick('log_scale', "Plot on a Log Scale?", master, ROW, COL)
   #COL+=1
   #tick('ThreeD', "Plot the Intensity in 3D", master, ROW, COL)
   #COL-=1
   #ROW+=1
   #Label(master, text="Choose Alt/Az for 3D plots (degrees)").grid(row= ROW, column=COL, sticky = W)
   #ROW+=1
   #enter_num('altitude', "Altitude", master, ROW, COL)
   #ROW+=1
   #enter_num('azimuth', "Azimuth", master, ROW, COL)
   #ROW+=1
   enter_num('pixels', "Number of Pixels (x y)", master, ROW, COL)
   ROW+=1
   enter_num('ave_dist', "Neighbouring Point Distance (nm)", master, ROW, COL)
   ROW+=1
   enter_num('z_scale','z-direction scaling', master, ROW, COL)
   ROW+=1
   tick('bound', "Upper / Lower Bounds?", master, ROW, COL)
   #g.dictionary_in['bound2']['font'] = "Times 11 underline"
   ROW+=1
   enter_num('minimum', "Minimum (ie. 3e-7)", master, ROW, COL)
   ROW+=1
   enter_num('maximum', "Maximum", master, ROW, COL)
   ROW+=1
   



   ### Sequence Options ###
   
   COL+=2
   ROW=0
   
   Button(master, text='Real Space', command=plot_points, font = "Times 16 bold").grid(row=ROW, column = COL, sticky=W, pady=2)
   ROW+=1
   
   #if g.dictionary['seq_hide'] == 0:
   #   int_button = Button(master, text="Calculate Intensity", command = calc_int, font = "Times 16 bold")
   #else:
   #   int_button = Button(master, text="Calculate Sequence", command = calc_seq, font = "Times 16 bold")
   int_button = Button(master, text="Calculate Intensity", command = calc_int, font = "Times 16 bold")
   int_button.grid(row=ROW, column = COL,sticky=W, pady=2)
   ROW+=1
   
   Button(master, text='Replot Intensity', command=view_intensity, font = "Times 16 bold").grid(row=ROW, column = COL, sticky=W, pady=2)
   
#   ROW+=1
#   Button(master, text='Plot a Ring', command=circ, font = "Times 16 bold").grid(row=ROW, column = COL,sticky=W, pady=2)
   
   ROW+=2
   Label(master, text="Pop-Up Windows:", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   ROW+=1

   Button(master, text='Ring Options', command=ring_options, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)   
   ROW+=1

   Button(master,text='Inter-Particle Scattering', command=interparticle_options, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)

   ROW+=1
   
   Button(master, text='Detector Options', command=detector_parameters, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)   
   ROW+=1
   
   Button(master, text='Sequence Options', command=sequence_parameters, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)   
   ROW+=1
   
   # COL += 2
   # ROW = 0
   # Label(master, text="Sequence Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   # COL+=1
   # if g.dictionary['seq_hide'] == 0:
   #    seq_button = Button(master, text="Edit Sequence", font = "Times 12 bold")
   # else:
   #    seq_button = Button(master, text="No Sequence", font = "Times 12 bold")
   # seq_button.bind("<Button-1>", hide_sequence)
   # seq_button.grid(row=ROW, column = COL, pady=2)
   # COL-=1
   # ROW+=1
   # enter_num('s_step', 'Number of Frames', master, ROW, COL)
   # ROW+=1
   # enter_str('s_var', "Which Variable?", master, ROW, COL)
   # ROW+=1
   # Button(master, text='Common Variables', command=show_sequence_variables).grid(row=ROW, column = COL, sticky=W, pady=2)
   #
   # ROW+=1
   # MODES_gauss = [('Linear Sequence', '0'),('Gaussian', '1'),]
   # radio('gauss',MODES_gauss, master, ROW, COL)
   # ROW+=1
   # enter_num('s_start', "Sequence Start", master, ROW, COL)
   # ROW+=1
   # enter_num('s_stop', "Sequence Stop", master, ROW, COL)
   # ROW+=1
   # enter_num('SD', 'Standard Deviation', master, ROW, COL)
   # ROW+=1



   ### Fitting Options ###

   Button(master, text='Fitting Options', command=select_fit_parameters, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)
   ROW+=1

   Button(master, text='More Output Options', command=output_options, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)   
   ROW+=1


   if g.debug:
      ROW+=1
      Button(master, text='Debug Window', command=open_debug, font = "Times 14 bold").grid(row=ROW, column = COL, sticky=W, pady=2)   
      ROW+=1

#   Label(master, text="Fitting Options", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
#   ROW += 1
#   enter_num('fit_file', "Experimental Data Filename", master, ROW, COL)
#   ROW += 1
#   enter_num('center', "Center of Beamstop (x y)", master, ROW, COL)
#   ROW += 1
#   enter_num('border', "Additional Cropping", master, ROW, COL)
#   ROW += 1
#   enter_num('mask_threshold', "Mask Threshhold", master, ROW, COL)
#   ROW += 1
#   COL += 1
#   Button(master, text="Plot Exp Data", command = plot_exp_data, font = "Times 16 bold").grid(row=ROW, column=COL, pady=2)
#   COL -= 1
#   ROW += 1
#   Label(master, text="Fit Parameters:").grid(row= ROW, column=COL, columnspan =2, sticky = W)
#   ROW += 1
#   enter_num('background', "Background Noise", master, ROW, COL)
#
#   ROW+=1      # These are checkboxes which, if unchecked, will hold fixed fit parameters.
#   tick("fit_radius_1", "Radius 1", master,ROW,COL)
#   COL+=1
#   tick("fit_radius_2", "Radius 2", master,ROW,COL)
#   COL-=1
#   ROW+=1
#   tick("fit_rho_1", "Rho 1", master,ROW,COL)
#   COL+=1
#   tick("fit_rho_2", "Rho 2", master,ROW,COL)
#   COL-=1
#   ROW+=1
#   tick("fit_z_dim", "Length", master,ROW,COL)
#   COL+=1
#   tick("fit_x_theta", "x rotation", master,ROW,COL)
#   COL-=1
#   ROW+=1
#   tick("fit_y_theta", "y rotation", master,ROW,COL)
#   COL+=1
#   tick("fit_z_theta", "z rotation", master,ROW,COL)
#   COL-=1
#   ROW+=1
#   tick("fit_background", "background", master,ROW,COL)
#   COL+=1
#   tick("fit_other", "unused", master,ROW,COL)
#   COL-=1
#   ROW += 1
#   ROW += 1
#   enter_num('max_iter', "Maximum Iterations (0=default)", master, ROW, COL)
#   ROW += 1
#   if g.debug:
#      enter_num('update_freq', "Update Interval", master, ROW, COL)   #TODO: debug
#      ROW += 1
#   enter_num('grid_compression', "Grid Compression (2, 5, or 10)", master, ROW, COL)
#   ROW += 1
#   tick('plot_fit_tick',"Plot Fit Results", master,ROW, COL)
#   COL += 1
#   tick('plot_residuals_tick',"Plot Fit Residuals", master,ROW, COL)
#   COL-=1
#   ROW += 1
#   Button(master, text="Plot Residuals", command = plot_residuals, font = "Times 16 bold").grid(row=ROW, column=COL, pady=2)
#   COL+= 1
#   Button(master, text="Fit Exp Data", command = perform_fit, font = "Times 16 bold").grid(row=ROW, column=COL, pady=2)
#   COL -= 1

   #select_fit_parameters()
   
   rename_parameters(0)
   root.mainloop()
   



