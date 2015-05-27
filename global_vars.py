#list of global variables


opencl_enabled = True
f2py_enabled = True
accelerate_points = True   #In case there are bugs and you want to disable only this portion of the acceleration.

quiet = False
verbose = False
debug = False

#default settings:
dictionary = {'advanced':1, 'altitude':45, 'analytic': 2, 'ave_dist': 0.6, 'azimuth':45, 'bound': 1, 'circ_delta':5, 'comments':'',
              'degrees': 1, 'energy_wavelength': 12, 'energy_wavelength_box': 0, 'gauss':0, 'log_scale': 1, 'maximum': 0.01, 'minimum': 1e-8,
              'num_plots': 1, 'pixels': 200, 'proportional_radius':0.5, 'QSize': 6,'Qz': 0, 'radius_1': 5.0, 'radius_2': 2.5, 'rho_1': 1.0, 'rho_2': -0.5,
              'save_img':1, 'save_name': 'save_name', 'scale': 1,'SD':1, 'seq_hide':0, 'shape': 2, 's_start': 0, 's_step': 2,
              's_stop': 1, 'subfolder':'subfolder', 's_var': 'x_theta', 'symmetric': 0,
              'theta_delta':20, 'ThreeD': 0, 'title': 'title', 'x_theta': 0,'y_theta': 0,'z_theta': 0,'z_dim': 100,'z_scale':1,
              'fit_file': 'fit_file', 'center': (0,0), 'border': 0, 'max_iter': 1000, 'update_freq': 0, 'plot_fit_tick': 1, 'plot_residuals_tick': 1, 'mask_threshold': 10, 'background': 2e-5, 'grid_compression': 5,
              'fit_radius_1': 1, 'fit_radius_2': 0, 'fit_rho_1': 1, 'fit_rho_2': 0, 'fit_z_dim': 1, 'fit_x_theta': 1, 'fit_y_theta': 1, 'fit_z_theta': 1, 'fit_background': 1, 'fit_other': 0
              }

if debug:
   print("Debug mode is on.")
   quiet = False
   verbose = True

def vprint(x,level=1):
   try:
      if verbose >= level:
         print(x)
   except:
      return
   return

def dprint(x,level=2):
   try:
      if verbose >= level:
         print(x)
   except:
      return
   return

def lprint(text,filename=0):
   '''Simple routine for logging text printed to the screen.'''
   print(text)
   if filename:
      with open(filename,'a') as log:
         log.write(text)

