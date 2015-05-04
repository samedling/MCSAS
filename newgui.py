#!/Users/scott/Library/Enthought/Canopy_64bit/User/bin/python
#!/usr/bin/python

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
#from scipy import ndimage  #possible smoothing of exp_data fore viewing

#Looks for fastmath.so to speed up intensity calculation.
try:
   import fastmath
   accelerated = True
   print "Accelerating using f2py."
except ImportError:
   accelerated = False
   print "Could not accelerate using f2py; if speedup is desired, run `make`."

quiet = False
verbose = False

#These are the default settings
dictionary = {'advanced':1, 'altitude':45, 'analytic': 2, 'ave_dist': 0.6, 'azimuth':45, 'bound': 1, 'circ_delta':5, 'comments':'',
              'degrees': 1, 'energy_wavelength': 12, 'energy_wavelength_box': 0, 'gauss':0, 'log_scale': 1, 'maximum': 0.01, 'minimum': 1e-8,
              'num_plots': 1, 'pixels': 200, 'proportional_radius':0.5, 'QSize': 6,'Qz': 0, 'radius_1': 5.0, 'radius_2': 2.5, 'rho_1': 1.0, 'rho_2': -0.5,
              'save_img':1, 'save_name': 'save_name', 'scale': 1,'SD':1, 'seq_hide':0, 'shape': 2, 's_start': 0, 's_step': 2,
              's_stop': 1, 'subfolder':'subfolder', 's_var': 'x_theta', 'symmetric': 0,
              'theta_delta':20, 'ThreeD': 0, 'title': 'title', 'x_theta': 0,'y_theta': 0,'z_theta': 0,'z_dim': 100,'z_scale':1,#}
              'fit_file': 'fit_file', 'center': (0,0), 'border': 0, 'max_iter': 1000, 'update_freq': 0, 'plot_fit_tick': 1, 'plot_residuals_tick': 1, 'mask_threshold': 10, 'background': 2e-5, 'grid_compression': 5,
              'fit_radius_1': 1, 'fit_radius_2': 0, 'fit_rho_1': 1, 'fit_rho_2': 0, 'fit_z_dim': 1, 'fit_x_theta': 1, 'fit_y_theta': 1, 'fit_z_theta': 1, 'fit_background': 1, 'fit_other': 0
              }
length_dictionary = len(dictionary)

#####            Importing data or using defaults              #############

#root_folder = os.path.dirname(sys.argv[0]) #Doesn't work when called from ipython.
root_folder = os.getcwd()
#print root_folder
#Check for write access?

try:
    d = pickle.load(open(root_folder+"/default.txt", 'rb'))
    if length_dictionary != len(d): #I check that it is the same length - This is needed if any new variables are added to dictionary
        a= 1/0
    dictionary = d
except:
    print "Previously used variables could not be loaded. \nUsing default settings instead."
    with open(root_folder+"/default.txt", 'wb') as f:
       pickle.dump(dictionary, f)

dictionary = {x:dictionary[x] for x in dictionary} #This contains the unaltered parameters
dictionary_in = {a:dictionary[x] for x in dictionary for a in [x, x+'2']} #This contains raw data from the GUI. The +'2' is so that i can control checkboxes
dictionary_SI = {x:dictionary[x] for x in dictionary} #this dictionary has the parameters after they have been converted to SI units

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
                        ["Double Cone",14]
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
   global dictionary, dictionary_SI
   dictionary_SI['x_dim'] = 2.*dictionary_SI['radius_1']
   dictionary_SI['y_dim'] = 2.*dictionary_SI['radius_1']
   
   

def get_numbers_from_gui():
    global dictionary, dictionary_in, dictionary_SI
    #Here I get all parameters from the GUI and put them into dictionary
    for x in dictionary:
       if x=='comments':
          dictionary[x] = dictionary_in[x].get(1.0,END).rstrip()
       elif x=='advanced' or x== 'seq_hide':
          dictionary[x] = dictionary[x]
       elif x=='shape':
          dictionary[x] = MC_num_and_name_dict[dictionary_in['shape'].get()]
       elif x=='analytic':
          dictionary[x] = Analytic_dict[dictionary_in['analytic'].get()]
       else:
          dictionary[x] = dictionary_in[x].get() 
    
    for x in dictionary:
        try:
            dictionary[x] = float(dictionary[x]) #I am turning the numbers from strings to floats.
            if dictionary[x]==int(dictionary[x]):#I am turning floats to integers if they are the same
               dictionary[x] = int(dictionary[x])
        except:
            None
    if not os.path.exists(root_folder+'/'+dictionary_SI['subfolder']):#making the subfolder, if it doesn't exist
       os.makedirs(root_folder+'/'+dictionary_SI['subfolder'])
       time.sleep(2) #Making the subfolder takes a few seconds, so we need to delay the program, otherwise it will try save things into the folder before it is made.

    dictionary_SI = {x: dictionary[x] for x in dictionary}
    dictionary_SI['path_to_subfolder'] = os.path.join(root_folder,dictionary['subfolder'],dictionary['save_name']) #This is for convienience

    #Converting to SI units.
    dictionary_SI["z_dim"] = dictionary["z_dim"]*10**-9
    dictionary_SI["ave_dist"] = dictionary["ave_dist"]*10**-9
    dictionary_SI["travel"] = dictionary_SI["ave_dist"]
    dictionary_SI["radius_1"] = dictionary["radius_1"]*10**-9
    xy_dim()#defining x_dim and y_dim - dependent of radius_1
    dictionary_SI["radius_2"] = dictionary["radius_2"]*10**-9
    dictionary_SI["QSize"] = dictionary["QSize"]*10**9

    dictionary_SI['num_plot_points'] = int(dictionary_SI['pixels']/2.)
    dictionary_SI['delta'] = 1. #number of pixels in width


    dictionary_SI['circ_delta'] = 1.4*dictionary_SI['circ_delta']
    dictionary_SI['pixel_radius'] = dictionary['proportional_radius']*dictionary_SI['pixels']/2.
    dictionary_SI['theta_delta'] = 6.283/dictionary['theta_delta']
    


    if dictionary_SI["energy_wavelength_box"] == 0: #Checkbox
       dictionary_SI["EHC"] = 2.*np.pi*(dictionary_SI["energy_wavelength"])*1.602176487*10**10/(6.62606896*2.99792458) #Energy to 2pi/lambda
    else:
       dictionary_SI["EHC"] = 2. * np.pi / dictionary_SI["energy_wavelength"]#lambda to 2pi/lambda
    if dictionary_SI["degrees"] == 1: #Conveting to radians
       dictionary_SI["x_theta"] = dictionary["x_theta"]*np.pi/180
       dictionary_SI["y_theta"] = dictionary["y_theta"]*np.pi/180
       dictionary_SI["z_theta"] = dictionary["z_theta"]*np.pi/180

    with open(root_folder+"/default.txt", 'wb') as f:
        pickle.dump(dictionary, f)#Saving the infomation from dictionary so it can be loaded later

    with open(dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
       pickle.dump(dictionary, f) #saving a copy in the subfolder for reference.


def load_functions(): #This loads the functions from the other files. It needs to be dynamic, hence I cannot use import.
   global dictionary, dictionary_SI
   execfile(root_folder+"/Monte_Carlo_Functions.py",globals())
   execfile(root_folder+"/Plotting_Functions.py", globals())
   execfile(root_folder+"/density_formula.py", globals())
   execfile(root_folder+"/analytic_formula.py", globals())

def change_units(number): #Used for sequences. A value is converted to SI units.
        global dictionary_SI
        for x in dictionary_SI:
           if x == dictionary_SI['s_var']:
              dictionary_SI[x] = number*10**-9
              if x == 'ave_dist':
                 dictionary_SI['travel'] = dictionary_SI[x]
              if x == 'travel':
                 dictionary_SI['ave_dist'] = dictionary_SI[x]
              if dictionary_SI['s_var'] != 'y_dim' and dictionary_SI['s_var'] != 'x_dim':
                 xy_dim()#This is here, mainly for the double slit - if you want to make the slits higher, you can. Most other functions are radially symmetric.
                               
              if x == "QSize":
                 dictionary_SI[x] = number*10**9
              if x == "energy_wavelength":
                  if dictionary_SI['energy_wavelength_box'] == 1:
                     dictionary_SI['EHC'] = 2.*np.pi*number*1.602176487*10**10/(6.62606896*2.99792458)
                  else:
                     dictionary['EHC'] = 2. * np.pi / number
              if x == "x_theta" or x == "y_theta" or x == "z_theta":
                  if dictionary_SI['degrees'] == 1:
                      dictionary_SI[x] = number*np.pi/180

def save_vars_to_file(extra): #here I save all the infomation into a text file that is easy to read. extra is a string of extra infomation that you might want to include.
   global dictionary_SI, dictionary
   get_numbers_from_gui()
   sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
   sim_info.write("\n\n\nDATE: "+time.strftime("%d")+time.strftime("%b")+time.strftime("%y")+"    TIME: "+time.strftime("%X")+"\n")
   sim_info.write(dictionary_SI['comments'])
   print dictionary_SI['path_to_subfolder']
   print dictionary['comments']
   print extra
   sim_info.write(extra+"\n")
   for x in sorted(dictionary, key=lambda x: x.lower()):
      if x != 'comments':
         sim_info.write(x+" : " +str(dictionary[x])+'\n')
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
    global dictionary, dictionary_SI, dictionary_in
    get_numbers_from_gui()
    save_vars_to_file("Plot Points")
    load_functions()
    if dictionary['seq_hide'] == 1:
       if dictionary['gauss']==0:
          current_value = dictionary['s_start']
       else:
          current_value = np.random.normal(loc = dictionary[dictionary['s_var']], scale = dictionary['SD'])
       change_units(current_value)
    Points_Plot(Points_For_Calculation(), 'points', 1)
    clear_mem()
    print "Program Finished"

def view_intensity(): #This allows you to view a premade intensity
    global dictionary, dictionary_SI, dictionary_in
    get_numbers_from_gui()
    load_functions()
    radial_intensity = pylab.loadtxt(dictionary_SI['path_to_subfolder']+"radial_intensity.csv", delimiter=",")
    radial_intensity_plot(radial_intensity, "radial", dictionary_SI['title'], 0)
    Intensity = pylab.loadtxt(dictionary_SI['path_to_subfolder']+"intensity.csv", delimiter=",")
    Intensity_plot(Intensity, "intensity", dictionary_SI['title'], 1)
    clear_mem()
    print "Program Finished"
    
def make_intensity(): #This makes an intensity
    global dictionary, dictionary_SI, dictionary_in, sim_info
    get_numbers_from_gui()
    save_vars_to_file("Monte Carlo Intensity")
    load_functions()
    Intensity = Average_Intensity()
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    if dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    print "Program Finished"
    
def sequence(): #This makes a sequence of intensities
    global dictionary, dictionary_SI, dictionary_in, sim_info
    get_numbers_from_gui()
    save_vars_to_file("Monte Carlo Sequence")
    load_functions()
    for frame_num in range(int(dictionary['s_step'])):
       sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
       sim_info.write("\nFrame " + str(frame_num+1) + " of " + str(int(dictionary['s_step'])))
       sim_info.close()

       print "\nmaking frame " + str(frame_num+1) + " of " + str(int(dictionary['s_step']))
       if dictionary['gauss']==0:
          try:
             current_value = (dictionary['s_stop']-dictionary['s_start'])*frame_num/(1.*dictionary['s_step']-1.)+1.*dictionary['s_start']
          except ZeroDivisionError:
             current_value = dictionary['s_start']
       else:
          current_value = np.random.normal(loc = dictionary[dictionary['s_var']], scale = dictionary['SD'])
       change_units(current_value)
       Intensity = Average_Intensity()
       save(Intensity, "intensity"+str(frame_num+1))
       radial_intensity = radial(Intensity)
       save(radial_intensity, "radial_intensity"+str(frame_num+1))
       try:
           cumulative += np.asarray(Intensity)
       except NameError:
           cumulative = np.asarray(Intensity)
       title = dictionary_SI['title']+" "+dictionary_SI['s_var']+'='+str(current_value)
       if dictionary_SI['save_img'] == 1:
          Intensity_plot(Intensity, "intensity" + str(frame_num+1), title, 0)
          radial_intensity_plot(radial_intensity, "radial_intensity" + str(frame_num+1), title, 0)
       clear_mem()
    Intensity = cumulative / dictionary_SI['s_step']
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    dictionary_SI['title'] = dictionary_SI['title']+" Averaged" + dictionary_SI['s_var']
    if dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    print "Program Finished"


def theory_plot(): #This plots an analytic model
   global dictionary_SI
   get_numbers_from_gui()
   save_vars_to_file("Analytic Intensity")
   load_functions()
   Intensity = theory_csv()
   save(Intensity, "intensity")
   radial_intensity = radial(Intensity)
   save(radial_intensity, "radial_intensity")
   if dictionary_SI['save_img'] == 1:
      view_intensity()
   clear_mem()
   print "Program Finished"



def theory_seq(): #This plots a sequence created with the analytic model
    global dictionary, dictionary_SI, dictionary_in, sim_info
    get_numbers_from_gui()
    save_vars_to_file("Analytic Sequence")
    load_functions()
    for frame_num in range(int(dictionary['s_step'])):
       sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
       sim_info.write("\nFrame " + str(frame_num+1) + " of " + str(int(dictionary['s_step'])))
       sim_info.close()

       print "\nmaking frame " + str(frame_num+1) + " of " + str(int(dictionary['s_step']))
       if dictionary['gauss']==0:
          try:
             current_value = (dictionary['s_stop']-dictionary['s_start'])*frame_num/(1.*dictionary['s_step']-1.)+1.*dictionary['s_start']
          except ZeroDivisionError:
             current_value = dictionary['s_start']
       else:
          current_value = np.random.normal(loc = dictionary[dictionary['s_var']], scale = dictionary['SD'])
       change_units(current_value)
       Intensity = theory_csv()
       save(Intensity, "intensity"+str(frame_num+1))
       radial_intensity = radial(Intensity)
       save(radial_intensity, "radial_intensity"+str(frame_num+1))
       try:
           cumulative += np.asarray(Intensity)
       except NameError:
           cumulative = np.asarray(Intensity)
       title = dictionary_SI['title']+" "+dictionary_SI['s_var']+'='+str(current_value)
       if dictionary_SI['save_img'] == 1:
          Intensity_plot(Intensity, "intensity" + str(frame_num+1), title, 0)
          radial_intensity_plot(radial_intensity, "radial_intensity" + str(frame_num+1), title, 0)
       clear_mem()
    Intensity = cumulative / dictionary_SI['s_step']
    save(Intensity, "intensity")
    radial_intensity = radial(Intensity)
    save(radial_intensity, "radial_intensity")
    dictionary_SI['title'] = dictionary_SI['title']+" Averaged" + dictionary_SI['s_var']
    if dictionary_SI['save_img'] == 1:
      view_intensity()
    clear_mem()
    print "Program Finished"



def circ(): #This plots a the angle at a fixed radius
   get_numbers_from_gui()
   load_functions()
   Intensity = np.asarray(pylab.loadtxt(dictionary_SI['path_to_subfolder']+"intensity.csv", delimiter=","))
   data = plotting_circle(Intensity)
   radial_intensity_plot(data, "theta"+str(dictionary['radius_2']), dictionary['title']+" "+str(dictionary['radius_2']), 0)
   angle_plot(data, "Angle"+str(dictionary['radius_2']), dictionary['title']+" "+str(dictionary['radius_2']), 1)
   print "finsihed"
   

def int_seq(): #This is the button, it runs a sequence or a single image depending on whether or not you can edit a sequence (For both analytic models and Monte Carlo Models)
   global dictionary_SI
   if dictionary['seq_hide'] == 0:
      if dictionary['shape'] ==0:
         theory_plot()
      else:
         make_intensity()
   else:
      if dictionary['shape']==0:
         theory_seq()
      else:
         sequence()







### Fitting Functions ###


class Fit_Parameters():
   '''Class to keep track of which parameters are being varied, based on shape and checkboxes.
      Contains functions to synchronize these parameters with global dictionaries.
      Also contains list of units for user-friendly output.'''
   def __init__(self):
      global dictionary,dictionary_SI
      shape=dictionary['shape']
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
         self.density_params=('radius_1','radius_2','rho_1')   #radius_1 only for scaling.
      elif shape in (8,11):   #Bubbles, Double Slit
         self.density_params=('radius_1','radius_2','rho_1')
      elif shape in (9,12,13):#Chopped Cylinder, N-Shaped Chopped Cone, Sine
         self.density_params=('radius_1','radius_2','rho_1','rho_2')  #z_dim intrinsic too
      elif shape == 10:
         print('Model not supported.')
      else:
         print('Unknown model. Assuming model uses all parameters.')
         self.density_params=('radius_1','radius_2','rho_1','rho_2')  #z_dim intrinsic too
      for name in (self.density_params+always):
         if dictionary['fit_'+name]:   #Looks for checkbox values.
            self.names.append(name)
      self.values=[dictionary_SI[var] for var in self.names]
      self.length=len(self.values)
      self.units=[]
      for i in range(self.length):     #Makes array of unit names correlated with values (for printing).
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

   def print_param(self,logfile=0):
      '''Prints parameters to the screen in a user-friendly format.'''
      global dictionary
      convert_from_SI()
      for i in range(self.length):
         lprint('{0} is {1:.4}{2}.'.format(self.names[i],float(dictionary[self.names[i]]),self.units[i]),logfile)
         #print('{0} is {1}{2}.'.format(self.names[i],self.values[i],self.units[i])) #Units always wrong?

   def get_param(self):
      '''Returns array usable by fitting routines.'''
      return self.values

   def set_param(self,parameters):
      '''Sets parameters from array from fitting routine.'''
      self.values = parameters

def load_exp_image(preview=False,enlarge_mask=1):
   '''Loads experimental data from file, cropping it and downsampling it if neccessary, and normalizes.  Also outputs the mask corresponding to the beamstop.'''
   #and 2D array containing a list of x,y points for the grid.'''
   global dictionary
   downsample=(dictionary['pixels'],dictionary['pixels'])
   center=[int(i) for i in dictionary['center'].split()]
   border=int(dictionary['border'])
   mask_threshold=dictionary['mask_threshold']
   filename=dictionary['fit_file']
   if preview:
      try:
         exp_data=np.array(Image.open(filename))
      except IOError:
         print('File {0} does not exist.'.format(filename))
         return
      normalize = 1.0/np.sum(exp_data)
      exp_data=exp_data*normalize
      #exp_data=ndimage.gaussian_filter(exp_data,sigma=3)
      return exp_data
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
      downsampled=cropped.resize(downsample,Image.BICUBIC)      #NEAREST,BILINEAR,BICUBIC,ANTIALIAS (worst to best; fastest to slowest; except ANTIALIAS does weird things sometimes)
      print("Resized to {0}.".format(downsample))
      exp_data=np.array(downsampled)
      padded=np.lib.pad(exp_data,((1,1),(1,1)),'edge')   #pads the array for enlarging mask
      mask=np.ones(np.product(exp_data.shape)).reshape(exp_data.shape)
      for i in range(mask.shape[0]):
         for j in range(mask.shape[1]):
            if exp_data[i,j] < mask_threshold:        #do after normalize?
               mask[i,j] = 0
            elif enlarge_mask and padded[i:i+3,j:j+3].min() < mask_threshold:        #do after normalize?
               mask[i,j] = 0  #could set to ~0.1 if want to decrease but not zero it.
      normalize = 1.0/np.sum(exp_data)
      exp_data=exp_data*normalize
      #img=Image.fromarray(exp_data)   #To go back to an image.
      return exp_data,mask

def plot_exp_data():#threshold=1e-7,zero_value=1e-7):
    '''Plots experimental data, using center and crop parameters if nonzero.'''
    global dictionary,dictionary_SI
    get_numbers_from_gui()
    load_functions()    #Needed for plotting routines.
    if dictionary['center'] == "0 0":
        image=load_exp_image(preview=True)
        print('Original image size is {0} x {1} pixels.'.format(image.shape[0],image.shape[1]))
        print('This takes a minute...')
        threshold=np.median(image)/10
        zero_value=threshold
        image[image<threshold]=zero_value
        Intensity_plot(image,"exp_data",'Before Cropping or Downsampling',1)
    else:
        cropped=load_exp_image()
        #cropped[cropped<threshold]=zero_value
        #Intensity_plot(cropped,"exp_data2",'After Cropping and Downsampling',1)
        #threshold=np.median(cropped[0])/10      #Remove this since upper/lower bounds is the same.
        #zero_value=threshold
        masked=cropped[0]*cropped[1]
        #masked[masked<threshold]=zero_value
        #for i in (2,5,10,20,50,80,90,95,98,99):
            #print('{0}th percentile value is {1:0.4}.'.format(i,np.percentile(masked,i)))
        print('Recommended value for background noise parameter is {0:0.4}.'.format(np.percentile(masked,25)))
        Intensity_plot(masked,"exp_data2",'After Cropping and Downsampling',1)
    return



def residuals(param,exp_data,mask=[],random_seed=2015):
   '''Returns residual array of difference between experimental data and data calculated from passed parameters.'''
   global dictionary_SI
   global parameters
   parameters.set_param(param)
   parameters.sync_dict()
   x=range(exp_data.shape[0])
   y=range(exp_data.shape[1])
   if not len(mask):
      mask = np.ones(exp_data.shape)
   err = np.zeros(np.product(exp_data.shape)).reshape(exp_data.shape)
   #load_functions()    #DO I NEED?  #Reintilizes functions with the new parameters.
   #calc_intensity = Average_Intensity() #might just take longer or might be necessary to accomodate randomness in Points_For_Calculation
   calc_intensity = Detector_Intensity(Points_For_Calculation(seed=random_seed),mask)  #like Average_Intensity() but just runs once and without time printouts and with same random_seed
   for i in x:
      for j in y:
         #err[i,j] = exp_data[i,j]-calc_intensity[i,j]
         if mask[i,j]:
            err[i,j] = mask[i,j]*(exp_data[i,j]-(calc_intensity[i,j]+dictionary_SI['background']))
         #err[i,j] = exp_data[i,j]-max(calc_intensity[i,j],dictionary_SI['background'])
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
   global dictionary_SI,parameters
   guess = parameters.get_param()
   fit_param = leastsq(residuals,guess,args=(exp_data,mask,int(time.time())),full_output=1,maxfev=update_freq)
   #parameters.set_param(fit_param[0])   #These two lines shouldn't really be needed.
   #parameters.sync_dict()
   with open(root_folder+"/default.txt", 'wb') as f:
      pickle.dump(dictionary, f)#Saving the infomation from dictionary so it can be loaded later
   with open(dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
      pickle.dump(dictionary, f) #saving a copy in the subfolder for reference.
   return fit_param

def perform_fit():  #Gets run when you press the Button.
   '''Loads experimental data from filename, fits the data using current dictionary as initial guesses, leaves final parameters in dictionary.'''
   global dictionary,dictionary_SI,parameters,quiet
   get_numbers_from_gui()
   load_functions()
   filename = dictionary['fit_file']
   max_iter = dictionary['max_iter']
   update_freq = dictionary['update_freq']
   plot_fit=dictionary['plot_fit_tick']
   plot_diff=dictionary['plot_residuals_tick']
   logfile=dictionary_SI['path_to_subfolder']+'fitlog.txt'   #dictionary['fitlog']
   grid_compression=dictionary['grid_compression']
   if not accelerated:
      print('Fortran acceleration is NOT enabled!')
      if grid_compression > 1:
         print('Grid compression does not work without Fortran.')
      print('This will probably take a REALLY LONG time.')
   initial_quiet=quiet
   quiet = True
   total_steps = 0
   if update_freq == 0:
      update_freq = max_iter
   exp_data,mask=load_exp_image()
   parameters=Fit_Parameters()  #Creates class of parameters with values and names.
   lprint('Starting values:',logfile)
   parameters.print_param(logfile)
   lprint('{0}: Starting fit...'.format(time.strftime("%X")),logfile)
   if grid_compression:
      mask = fast_mask(exp_data,mask,grid_compression)
   total_steps = 0
   while total_steps < max_iter or max_iter == 0:
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
   quiet=initial_quiet
   #need to refresh dictionary_SI?
   save(diff,"_fit_residuals")
   #with open(root_folder+"/default.txt", 'wb') as f:
   #   pickle.dump(dictionary, f)#Saving the infomation from dictionary so it can be loaded later
   #with open(dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
   #   pickle.dump(dictionary, f) #saving a copy in the subfolder for reference.
   if plot_fit and plot_diff:
      fit_results=Average_Intensity()
      save(fit_results,"_fit")
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
   '''Speeds calculation by adding points to mask.  Speedup can be (2,5,10).'''
   if speedup < 2:
      return mask
   elif speedup == 2:
      mod=2
      percentile=80
      pad=1
   elif speedup == 5:
      mod=3
      percentile=94
      pad=1
   elif speedup == 10:
      mod=5
      percentile=95
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
   print('Of {0} pixels, {1} are masked by beamstop and {2} are being skipped for speed.'.format(total,int(starting),int(final-starting)))
   print('Estimated speedup: {0:.3}x.'.format(total/(total-final)))
   return mask

def convert_from_SI():
    '''Copies parameters from SI dictionary back to regular dictionary.'''
    global dictionary,dictionary_SI
    dictionary["radius_1"] = dictionary_SI["radius_1"]*10**9
    dictionary["radius_2"] = dictionary_SI["radius_2"]*10**9
    dictionary["z_dim"] = dictionary_SI["z_dim"]*10**9
    if dictionary["degrees"] == 1: #Converting from radians
       dictionary["x_theta"] = dictionary_SI["x_theta"]*180/np.pi
       dictionary["y_theta"] = dictionary_SI["y_theta"]*180/np.pi
       dictionary["z_theta"] = dictionary_SI["z_theta"]*180/np.pi
    dictionary["rho_1"] = dictionary_SI["rho_1"]
    dictionary["rho_2"] = dictionary_SI["rho_2"]

def plot_residuals():
   '''Loads exp data, calculates intensity, and plots the difference [as well as 2 original plots].'''
   global dictionary
   plot_all=dictionary['plot_fit_tick']
   get_numbers_from_gui()
   load_functions()
   filename = dictionary['fit_file']
   if not accelerated:
      print('Fortran acceleration is NOT enabled!')
   print('{0}: Starting calculation...'.format(time.strftime("%X")))
   exp_data,mask=load_exp_image()
   calc_intensity=Average_Intensity()
   save(calc_intensity,"_calc")     #wrong suffix!!
   #err = np.zeros(np.product(exp_data.shape)).reshape(exp_data.shape)
   err = np.zeros(exp_data.shape)
   for i in range(exp_data.shape[0]):
      for j in range(exp_data.shape[1]):
         #err[i,j] = exp_data[i,j]-calc_intensity[i,j]
         if mask[i,j]:
            err[i,j] = mask[i,j]*(exp_data[i,j]-(calc_intensity[i,j]+dictionary['background']))
   save(err,"_guess_residuals")
   plot_residuals=np.abs(err)
   print('{0}: Total error = {1:.4}; sum of squares = {2:.4}'.format(time.strftime("%X"),plot_residuals.sum(),np.square(err).sum()))
   if plot_all:
      #view_fit(exp_data*mask,calc_intensity,plot_residuals)
      Fit_plot(exp_data*mask,calc_intensity,plot_residuals)
   else:
      print('Plotting difference.')
      Intensity_plot(plot_residuals,"residuals",'Difference Plot',1)

#Unused:
def view_fit(exp_data,fit_results,fit_residuals):
   '''Copied from view_intensity() with minor changes to plot all relevant fitting plots.'''
   global dictionary,dictionary_SI
   plot_fit = dictionary['plot_fit_tick']
   plot_residuals = dictionary['plot_residuals_tick'] #local name
   #get_numbers_from_gui()
   #load_functions()
   if plot_fit:
      print('Plotting exp data and calculated results.')
      threshold=np.median(exp_data)/10
      zero_value=threshold
      exp_data[exp_data<threshold]=zero_value
      Intensity_plot(exp_data,"exp_data",'Experimental Data',1)
      Intensity_plot(fit_results,"fit",'Calculated Data',1)
   if plot_residuals:
      print('Plotting difference.')
      Intensity_plot(fit_residuals,"residuals",'Difference Plot',1)
   clear_mem()
   print("Program Finished.")




def xy_grid(exp_data,percentile=95,mode=2,pad=1):
   '''Returns xy grid of points which are either above a threshold or every mode^2 datapoint.  Available modes are (2,3,5,10), corresponding to a decrease of (4,9,25,100) of points below threshold percentile.'''
   all_x=range(exp_data.shape[0])
   all_y=range(exp_data.shape[1])
   threshold=np.percentile(exp_data,percentile)
   x,y=[],[]
   if pad:
      padded=np.lib.pad(exp_data,((1,1),(1,1)),'edge')   #pads the array
      for i in all_x:
         for j in all_y:
            if padded[i:i+3,j:j+3].min() > threshold:
               x.append(i)
               y.append(j)
            elif mode == 10 and i%10==5 and j%10==5:  #factor <100
               x.append(i)
               y.append(j)
            elif mode == 5 and i%5==2 and j%5==2:     #factor <25
               x.append(i)
               y.append(j)
            elif mode == 3 and i%3==1 and j%3==1:     #factor <9
               x.append(i)
               y.append(j)
            elif mode == 2 and i%2==0 and j%2==0:     #factor <4
               x.append(i)
               y.append(j)
   else:
      for i in all_x:
         for j in all_y:
            if exp_data[i,j] > threshold:
               x.append(i)
               y.append(j)
            elif mode == 10 and i%10==5 and j%10==5:  #factor <100
               x.append(i)
               y.append(j)
            elif mode == 5 and i%5==2 and j%5==2:     #factor <25
               x.append(i)
               y.append(j)
            elif mode == 3 and i%3==1 and j%3==1:     #factor <9
               x.append(i)
               y.append(j)
            elif mode == 2 and i%2==0 and j%2==0:     #factor <4
               x.append(i)
               y.append(j)
   return x,y
### End Fitting Functions ###





################          MAKING THE GUI            #################
advanced = ['log_scale2', 'bound2', 'minimum', 'ave_dist', 'scale2',
              'ThreeD2', 'pixels', 'maximum', 'altitude', 'azimuth', 'z_scale']
#All these things will be hiden when advanced/simple options button is pressed

#This refers to 
def hide(event):
   if dictionary['advanced'] == 1:
      for x in advanced:
         dictionary_in[x].config(state=DISABLED)
      dictionary['advanced']=0
      advbutton["text"] = "Advanced Options"
   else:
      for x in advanced:
         dictionary_in[x].config(state=NORMAL)
      dictionary['advanced']=1
      advbutton["text"] = "Simple Options"


#These are all hidden when Edit/No sequence is pressed.
seq_options = ['s_step', 's_stop', 's_start', 's_var', 'gauss0', 'gauss1', 'SD']
def hide_sequence(event):
   if dictionary['seq_hide'] == 1:
      for x in seq_options:
         dictionary_in[x].config(state=DISABLED)
      dictionary['seq_hide']=0
      seq_button["text"] = "Edit Sequence"
      int_button['text'] = 'Calculate Intensity'
   else:
      for x in seq_options:
         dictionary_in[x].config(state=NORMAL)
      dictionary['seq_hide']=1
      seq_button["text"] = "No Sequence"
      int_button['text'] = 'Calculate Sequence'

#If you want a number box and a label, use this
def enter_num(variable_name, label, ROW, COL):
    global dictionary, dictionary_in
    Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
    dictionary_in[variable_name] = StringVar()
    dictionary_in[variable_name].set(dictionary[variable_name])
    dictionary_in[variable_name] = Entry(master, textvariable = dictionary_in[variable_name])
    if variable_name in advanced and dictionary['advanced'] == 0:
       dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and dictionary['seq_hide'] == 0:
       dictionary_in[variable_name].config(state = DISABLED)
    dictionary_in[variable_name].grid(row= ROW, column = COL+1)

#If you want a number box below the label
def enter_vert_num(variable_name, label, ROW, COL):
    global dictionary, dictionary_in
    Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
    dictionary_in[variable_name] = StringVar()
    dictionary_in[variable_name].set(dictionary[variable_name])
    dictionary_in[variable_name] = Entry(master, textvariable = dictionary_in[variable_name])
    if variable_name in advanced and dictionary['advanced'] == 0:
       dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and dictionary['seq_hide'] == 0:
       dictionary_in[variable_name].config(state = DISABLED)
    dictionary_in[variable_name].grid(row= ROW+1, column = COL)

#If you want a tick box
def tick(variable_name, label, ROW, COL):
    global dictionary, dictionary_in
    dictionary_in[variable_name] = IntVar()
    dictionary_in[variable_name].set(int(dictionary[variable_name]))
    dictionary_in[variable_name+'2'] = Checkbutton(master, text=label, variable=dictionary_in[variable_name])
    if variable_name+'2' in advanced and dictionary['advanced'] == 0:
       dictionary_in[variable_name+'2'].config(state = DISABLED)

    if variable_name+'2' in seq_options and dictionary['seq_hide'] == 0:
       dictionary_in[variable_name+'2'].config(state = DISABLED)
    dictionary_in[variable_name+'2'].grid(row=ROW, column = COL, sticky=W)

#If you want a string entered
def enter_str(variable_name, label, ROW, COL):
    global dictionary, dictionary_in
    Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
    dictionary_in[variable_name] = Entry(master)
    dictionary_in[variable_name].insert(0, dictionary[variable_name])
    if variable_name in advanced and dictionary['advanced'] == 0:
       dictionary_in[variable_name].config(state = DISABLED)
    if variable_name in seq_options and dictionary['seq_hide'] == 0:
       dictionary_in[variable_name].config(state = DISABLED)
    dictionary_in[variable_name].grid(row= ROW, column = COL + 1)

def enter_text(variable_name, label, WIDTH, HEIGHT, ROW, COL):#For a large textbox
   global dictionary, dictionary_in
   Label(master, text=label).grid(row= ROW, column=COL, sticky = W)
   dictionary_in[variable_name] = Text(master, height = HEIGHT, width = WIDTH)

   dictionary_in[variable_name].insert(1.0, dictionary[variable_name])
   dictionary_in[variable_name].grid(row=ROW+1, column = COL, rowspan = 2, sticky = W)

def radio(variable_name, MODES, ROW, COL): #Radiobutton
   global dictionary, dictionary_in
   dictionary_in[variable_name] = StringVar()
   dictionary_in[variable_name].set(dictionary[variable_name])
   for name, mode in MODES:
      dictionary_in[variable_name+mode] = Radiobutton(master, text=name, variable=dictionary_in[variable_name], value=mode)
      dictionary_in[variable_name+mode].grid(row=ROW, column = COL, sticky = W)
      if variable_name+mode in advanced and dictionary['advanced'] == 0:
          dictionary_in[variable_name+mode].config(state = DISABLED)
      if variable_name+mode in seq_options and dictionary['seq_hide'] == 0:
          dictionary_in[variable_name+mode].config(state = DISABLED)
      COL+=1




###########THE GUI STARTS HERE#####################

if __name__ == "__main__":
   master = Tk()
   master.title("Monte Carlo Small Angle Scattering (v 0.2.0), By Max Proft")

   ### Model Type ###

   ROW = 0
   COL = 0
   Label(master, text = "Model Type", font = "Times 16 bold").grid(row = ROW, column = COL, sticky = W)
   ROW+=1
   
   Label(master, text = "Choose a Monte Carlo Model").grid(row = ROW, column = COL, sticky = W)
   dictionary_in['shape'] = StringVar(master)
   dictionary_in['shape'].set(MC_num_and_name[dictionary['shape']][0])
   OptionMenu(master, dictionary_in['shape'], *MC_num_and_name[:,0]).grid(row = ROW, column = COL+1)
   
   ROW+=1
   tick("symmetric", "Radial Symmetry", ROW,COL)
   ROW+=1
   #COL+=1
   tick('Qz',"Small Angle Approx. (Qz=0)", ROW, COL)
   #COL-=1
   
   ROW+=1
   Label(master, text = "Choose an Analytic Model").grid(row = ROW, column = COL, sticky = W)
   dictionary_in['analytic'] = StringVar(master)
   dictionary_in['analytic'].set(Analytic_options[dictionary['analytic']][0])
   OptionMenu(master, dictionary_in['analytic'], *Analytic_options[:,0]).grid(row = ROW, column = COL+1)
   
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
   if dictionary['advanced'] == 0:
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
   enter_num('pixels', "Number of Pixels", ROW, COL)
   ROW+=1
   enter_num('ave_dist', "Neighbouring Point Distance (nm)", ROW, COL)
   ROW+=1
   enter_num('z_scale','z-direction scaling of\nneighbouring point distance', ROW, COL)
   ROW+=1
   tick('bound', "Upper and Lower Bounds?", ROW, COL)
   dictionary_in['bound2']['font'] = "Times 11 underline"
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
   if dictionary['seq_hide'] == 0:
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
   if dictionary['seq_hide'] == 0:
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
   enter_num('update_freq', "Update Interval", ROW, COL)
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
   



