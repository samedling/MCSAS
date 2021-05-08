import numpy as np
import global_vars as g
import copy

from scipy.optimize import curve_fit


class Fit_Parameters():
   '''Class to keep track of which parameters are being varied, based on shape and checkboxes.
      Contains functions to synchronize these parameters with global dictionaries.
      Also contains list of units for user-friendly output.'''
   def __init__(self):
      shape=g.dictionary['shape']
      
#      if shape in (1,2,6):    #Sphere, Cylinder, or Hex Prism
#         self.density_params=('radius_1','rho_1')
#      elif shape == 3:        #Core Shell
#         self.density_params=('radius_1','radius_2','rho_1','rho_2')
#      elif shape == 4:        #Gaussian
#         self.density_params=('radius_1','radius_2')     #radius_1 only for scaling.
#      elif shape in (5,14):   #Chopped Cone, Double Done
#         self.density_params=('radius_1','radius_2','rho_1') #z_dim intrinsic too
#      elif shape == 7:        #Rect. Prism
#         self.density_params=('radius_1','radius_2','rho_1')
#      elif shape in (8,11,15):   #Bubbles, Double Slit, Eliptical Cylinder
#         self.density_params=('radius_1','radius_2','rho_1')
#      elif shape in (9,12,13):#Chopped Cylinder, N-Shaped Chopped Cone, Sine
#         self.density_params=('radius_1','radius_2','rho_1','rho_2')  #z_dim intrinsic too
#      elif shape == 10:
#         print('Model not supported.')
#      else:
#         print('Unknown model. Assuming model uses all parameters.')
#         self.density_params=('radius_1','radius_2','rho_1','rho_2')  #z_dim intrinsic too
#      always=('x_theta','y_theta','z_theta','background','z_dim')
      
      #self.density_params=[g.var_list[i] if g.model_parameters[shape][2][i] for i in range(len(g.var_list))]
      self.density_params=[]
      for i in range(len(g.var_list)):
         if g.model_parameters[shape][2][i] != 'unused':
            self.density_params.append(g.var_list[i])
      always=['x_theta','y_theta','z_theta','background']

      self.names=[]
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
         elif self.names[i] in ('z_dim','length_2'):
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
         qsize_scale = (max(img.size)+0.)/max(r-l,t-b)
         g.dictionary_SI['QSize'] *= qsize_scale
         print('Detector Q Range has been increased to {0:4.4} due to cropping.'.format(g.dictionary['QSize']*qsize_scale))
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
      g.dprint("Normalizing to {0} so will be normalized to 1 when background of {1} is added.".format(norm_to,g.dictionary_SI['background']))
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
   #gprint("Normalized so total value is {0} and lowest value is {1}.".format(np.sum(normalized),np.min(normalized)))
   if g.verbose:
      g.dprint("Normalized so total value is {0} and lowest value is {1}.".format(np.sum(normalized*mask),np.min(normalized)))
   return normalized

def plot_exp_data():#threshold=1e-7,zero_value=1e-7):
    '''Plots experimental data, using center and crop parameters if nonzero.  Also shows the effect of grid compression if enabled.'''
    get_numbers_from_gui()
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
   global total_steps
   update_freq = g.dictionary['update_freq']
   parameters.set_param(param)
   parameters.sync_dict()
   x=range(exp_data.shape[0])
   y=range(exp_data.shape[1])
   if not len(mask):
      mask = np.ones(exp_data.shape)

   if g.dictionary['s_step'] > 1:    #sequence fit
      convert_from_SI()
      orig = g.dictionary[g.dictionary['s_var']]
      for i in range(g.dictionary['s_step']):
         #g.dictionary_SI[g.dictionary['s_var']] = np.random.normal(loc=g.dictionary_SI[g.dictionary['s_var']],scale=g.dictionary['SD'])
         if g.dictionary['gauss']:
            np.random.seed([random_seed+i])
            current_value = np.random.normal(loc=orig,scale=g.dictionary['SD'])
         else:
            start = g.dictionary['s_start']
            stop = g.dictionary['s_stop']
            current_value = (orig*2/(stop+start))*(start+((stop-start)*i/((g.dictionary['s_step']-1))))  #percent mode
            #current_value = (orig-(stop+start)/2)+(start+((stop-start)*i/((g.dictionary['s_step']-1))))  #difference mode
         change_units(current_value)
         g.vprint('{0} is {1} (SI: {2})'.format(g.dictionary['s_var'],current_value,g.dictionary_SI[g.dictionary['s_var']]))
         one_intensity = normalize(Calculate_Intensity(Points_For_Calculation(seed=random_seed),mask),mask,True)
         try:
            calc_intensity += one_intensity
         except NameError:
            calc_intensity = one_intensity
         #parameters.set_param(param)
         parameters.sync_dict()
   else:
      calc_intensity = normalize(Calculate_Intensity(Points_For_Calculation(seed=random_seed),mask),mask,True)

   weighted = True
   if weighted:
      #todo: if masked, actually set weight = 0?
      #sigma = sqrt(n) or sqrt(n)/n??!!
      g.dprint('Using weighted least squares fitting.')
      background=np.percentile(exp_data,20)
      sigma = np.sqrt(exp_data+background)               #sqrt(n)+background error bars
      #background=np.sqrt(np.percentile(exp_data,20))
      #sigma = np.maximum(np.sqrt(exp_data),background)   #sqrt(n) error bars (with 0 correction)
      #err = mask*(exp_data - calc_intensity)/sigma
      #err = np.log(mask*(exp_data - calc_intensity)/sigma)
   else:
      sigma = 1

   if g.dictionary['log_scale']:
      #TODO: doesn't fit properly
      #err = mask*np.log((exp_data / calc_intensity))/sigma
      #err[np.isnan(err)] = 0
      diff = exp_data / calc_intensity
      err = np.asarray([[ mask[i,j]*np.log(diff[i,j]) if (diff[i,j] > 0) else 0 for i in range(diff.shape[0]) ] for j in range(diff.shape[1])]) / sigma
   else:
      err = mask*(exp_data - calc_intensity)/sigma

   try:
      print('{0}: Step {3} total error = {1:.4}; sum of squares = {2:.4}'.format(time.strftime("%X"),np.abs(err).sum(),np.square(err).sum(),total_steps))
      total_steps +=1
      if update_freq and not total_steps%update_freq:
         logfile=g.dictionary_SI['path_to_subfolder']+'fitlog.txt'   #g.dictionary['fitlog']
         lprint('{0}: On function call {1}...'.format(time.strftime("%X"),total_steps),logfile,quiet=True)
         lprint('Current parameter values are:',logfile)
         parameters.print_param(logfile)
         if g.dictionary['log_scale']:
            if weighted:
               print('Fit type: weighted logarithmic.')
            else:
               print('Fit type: logarithmic.')
         else:
            if weighted:
               print('Fit type: weighted linear.')
            else:
               print('Fit type: linear.')
   except NameError:
      print('{0}: Total error = {1:.4}; sum of squares = {2:.4}'.format(time.strftime("%X"),np.abs(err).sum(),np.square(err).sum()))
   return np.ravel(err)     #flattens err since leastsq only takes 1D array

def lprint(text,filename=0,quiet=0):
   '''Simple routine for logging text printed to the screen.'''
   if not quiet:
      print(text)
   if filename:
      with open(filename,'a') as log:
         log.write(text+"\n")

def fit_step(exp_data,mask=[],max_iter=0):
   '''Runs a small number of iterations of fitting exp_data.'''
   global parameters
   guess = parameters.get_param()
   #norm_exp_data = exp_data/np.sum(exp_data*mask) #data normalized taking mask into account.
   exp_data = normalize(exp_data,mask)
   random_seed = int(time.time())

   fit_method = 'leastsq'
   #fit_method = 'weighted'
   if fit_method == 'weighted':  #TODO: doesn't work; for now use manual weighting in residuals.
      g.dprint('Using weighted least squares fitting.')
      fit_func = lambda params: residuals(params,exp_data,mask,random_seed)
      fit_sigma = np.ravel(exp_data)   #TODO: CHECK!!
      xdata=0
      ydata=0
      fit_param = curve_fit(fit_func,xdata,ydata,p0=guess,sigma=fit_sigma,absolute_sigma=True)
      #absolute_sigma=True means sigma are 1SD errors; otherwise, sigma are relative weights
   else:
      g.dprint('Using least squares fitting.')
      fit_param = leastsq(residuals,guess,args=(exp_data,mask,random_seed),full_output=1,maxfev=max_iter)

   with open(root_folder+"/default.txt", 'wb') as f:
      pickle.dump(g.dictionary, f)#Saving the infomation from g.dictionary so it can be loaded later
   with open(g.dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
      pickle.dump(g.dictionary, f) #saving a copy in the subfolder for reference.
   return fit_param

def perform_fit():  #Gets run when you press the Button.
   '''Loads experimental data from filename, fits the data using current g.dictionary as initial guesses, leaves final parameters in g.dictionary.'''
   global parameters
   get_numbers_from_gui()
   if g.dictionary['center'] == '0 0':
      print('ERROR: You must specify the centre of the beamstop before fitting.')
      print('Click Plot Exp. Data, hover the mouse over the centre, and record the x,y values.')
      return
   filename = g.dictionary['fit_file']
   max_iter = g.dictionary['max_iter']
   plot_fit=g.dictionary['plot_fit_tick']
   plot_diff=g.dictionary['plot_residuals_tick']
   logfile=g.dictionary_SI['path_to_subfolder']+'fitlog.txt'   #g.dictionary['fitlog']
   grid_compression=g.dictionary['grid_compression']
   #prefit_jumps = 5  #TODO: make into new user dictionary variable
   if max_iter:
      prefit_jumps = int(np.sqrt(max_iter))           #7-31 for 50-1000
      #prefit_jumps = int(2*np.power(max_iter,0.33))  #7-19 for 50-1000
      #prefit_jumps = int(3*np.power(max_iter,0.25))  #7-16 for 50-1000
   else:
      prefit_jumps=32
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
   global total_steps
   total_steps = 0
   exp_data,mask=load_exp_image()
   parameters=Fit_Parameters()  #Creates object of parameters with values and names.
   lprint('Starting values:',logfile)
   parameters.print_param(logfile)
   lprint('{0}: Starting fit...'.format(time.strftime("%X")),logfile)
   if grid_compression > 1:
      mask = fast_mask(exp_data,mask,grid_compression)
   if prefit_jumps:
      g.dprint("Starting {0} prefit random jumps.".format(prefit_jumps))
      jump_fit(exp_data,mask,iterations=prefit_jumps,jump_size=0.10)    #jump_size = 0.02 means +/- 1% change in each parameter
      g.dprint("Completed {0} prefit random jumps.".format(prefit_jumps))
      if g.debug:
         lprint('Prefit Values:',logfile)
         parameters.print_param(logfile)
   fit_param = fit_step(exp_data,mask,max_iter)
   #total_steps+=fit_param[2]['nfev']
   if fit_param[2]['nfev'] < max_iter:      #Checks if fit is completed.
      lprint('{0}: Converged after {1} function calls.'.format(time.strftime("%X"),total_steps),logfile)
      lprint('Final parameter values are:',logfile)
      parameters.print_param(logfile)
   else:
      lprint('{0}: Fit did not converge in {1} steps.'.format(time.strftime("%X"),total_steps),logfile)
      lprint('Parameter values are:',logfile)
      parameters.print_param(logfile)
   diff=fit_param[2]['fvec'].reshape(exp_data.shape)
   g.quiet=initial_quiet
   save(diff,"_fit_residuals")
   #with open(root_folder+"/default.txt", 'wb') as f:
   #   pickle.dump(g.dictionary, f)#Saving the infomation from g.dictionary so it can be loaded later
   #with open(g.dictionary_SI['path_to_subfolder']+"default.txt", 'wb') as f:
   #   pickle.dump(g.dictionary, f) #saving a copy in the subfolder for reference.
   if plot_fit and plot_diff:
      fit_results=normalize(Average_Intensity(),mask=mask,background=True)
      #fit_results=normalize(Calculate_Intensity(Points_For_Calculation()),background=True)
      save(fit_results,"_fit")
      diff = abs(exp_data - fit_results)*mask
      #diff = abs(exp_data - (fit_results + g.dictionary_SI['background']))
      Fit_plot(exp_data*mask,fit_results,diff)
   elif plot_fit:
      #fit_results=Average_Intensity()
      fit_results=normalize(Average_Intensity(),mask=mask,background=True)
      save(fit_results,"_fit")
      print('Plotting intensity.')
      #Intensity_plot(fit_results,"residuals",'Difference Plot',1)
      multiplot((exp_data*mask,fit_results),titles=('Exp. Data','Calc. Intensity'))
   elif plot_diff:
      print('Plotting difference.')
      Intensity_plot(diff,"residuals",'Difference Plot',1)


def jump_fit(exp_data,mask=[],vary_all=True,jump_size=0.02,init_residuals=0,seed=0,iterations=1):
   '''Randomly varies parameters (all, or a random one) if it lowers residuals.'''
   global parameters
   if not seed:
      seed = int(time.time())
   original_param = parameters.get_param()
   if not init_residuals:
      init_residuals = np.square(residuals(original_param,exp_data,mask,seed)).sum()
   np.random.seed([int(time.time())])
   vary_all = iterations%2 #TODO: decide if this is good!
   if vary_all:
      test_param = original_param*(1+(np.random.random(len(original_param))-0.5)*jump_size) #jump_size=0.02 results in up to +/-1% change in each parameter
   else: #vary only one, randomly, but by more
      test_param = copy.copy(original_param)
      n = np.random.randint(0,len(original_param))
      test_param[n] = original_param[n]*(1+(np.random.random()-0.5)*5*jump_size) #jump_size=0.02 results in up to +/-10% change in parameter
   test_residuals = np.square(residuals(test_param,exp_data,mask,seed)).sum()
   #np.random.seed([seed])
   #if test_residuals > init_residuals and np.random.random() > 0.1:  #10% chance of making jump anyway
   if test_residuals > init_residuals:
      parameters.set_param(original_param)
   else:
      g.dprint('New sum of squares is smaller; accepting new parameters.')
      init_residuals = test_residuals
   if iterations > 1:
      jump_fit(exp_data,mask,vary_all,jump_size,init_residuals,seed,iterations-1)


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
   get_numbers_from_gui()
   plot_fit=g.dictionary['plot_fit_tick']
   plot_diff=g.dictionary['plot_residuals_tick']
   filename = g.dictionary['fit_file']
   if not g.f2py_enabled and not g.opencl_enabled:
      print('Acceleration is NOT enabled!')
   print('{0}: Starting calculation...'.format(time.strftime("%X")))
   exp_data,mask=load_exp_image()
   if g.dictionary['grid_compression'] > 1:
      fast_mask(exp_data,mask,g.dictionary['grid_compression'])

   if g.dictionary['s_step'] > 1:    #sequence fit
      orig = g.dictionary[g.dictionary['s_var']]
      for i in range(g.dictionary['s_step']):
         #g.dictionary_SI[g.dictionary['s_var']] = np.random.normal(loc=g.dictionary_SI[g.dictionary['s_var']],scale=g.dictionary['SD'])
         if g.dictionary['gauss']:
            current_value = np.random.normal(loc=orig,scale=g.dictionary['SD'])
         else:
            start = g.dictionary['s_start']
            stop = g.dictionary['s_stop']
            current_value = (orig*2/(stop+start))*(start+((stop-start)*i/((g.dictionary['s_step']-1))))  #percent mode
            #current_value = (orig-(stop+start)/2)+(start+((stop-start)*i/((g.dictionary['s_step']-1))))  #difference mode
         change_units(current_value)
         g.vprint('{0} is {1} (SI: {2})'.format(g.dictionary['s_var'],current_value,g.dictionary_SI[g.dictionary['s_var']]))
         one_intensity = normalize(Calculate_Intensity(Points_For_Calculation(),mask),mask,True)
         try:
            calc_intensity += one_intensity
         except NameError:
            calc_intensity = one_intensity
   else:
      #calc_intensity = normalize(Calculate_Intensity(Points_For_Calculation(seed=random_seed),mask),mask,True)
      calc_intensity = normalize(Average_Intensity(mask),mask,True)  #Averages several and adds background.

   #calc_intensity = normalize(Average_Intensity(mask),mask,True)  #Averages several and adds background.

   #calc_intensity=Average_Intensity(mask)   #Old; probably doesn't add background
   #calc_intensity = normalize(Calculate_Intensity(Points_For_Calculation(),mask),mask,True) #Does not average.

   save(calc_intensity,"_calc")     #wrong suffix!!
   err = calc_intensity - exp_data
   ##calc_intensity = normalize(Calculate_Intensity(Points_For_Calculation(seed=random_seed),mask),mask,True)
   #err = mask*(exp_data - (calc_intensity + g.dictionary_SI['background']))
   #calc_intensity += g.dictionary_SI'background'] - exp_data  #todo: might be faster
   #calc_intensity *= mask                                     #todo: might be faster
   save(err,"_guess_residuals")
   plot_residuals=np.abs(err)*mask
   print('{0}: Total error = {1:.4}; sum of squares = {2:.4}'.format(time.strftime("%X"),plot_residuals.sum(),np.square(err).sum()))
   if plot_fit and plot_diff:
      Fit_plot(exp_data*mask,calc_intensity,plot_residuals)
   elif plot_fit:
      print('Plotting intensity.')
      multiplot((exp_data*mask,calc_intensity),titles=('Exp. Data','Calc. Intensity'))
   elif plot_diff:
      print('Plotting difference.')
      Intensity_plot(plot_residuals,"residuals",'Difference Plot',1)



def select_fit_parameters():
   global fit_window,master
   #fit_window = Tk()
   #fit_window.title("Select Fit Parameters")
   fit_window = Toplevel(master)
   fit_window.title("Select Fit Parameters")
   new_frame = Frame(fit_window)
   ROW = 0
   COL = 0
   #Label(fit_window, text = "Fit Parameters", font = "Times 16 bold").grid(row=ROW,column=COL,sticky=W)
   #COL+=1
   #Button(fit_window, text = "Save", command=save_fitparam, font= "Times 16 bold").grid(row=ROW,column=COL,pady=4)
   #COL-=1
   #ROW+=1
   

   ### Fitting Options ###

   Label(fit_window, text="Experimental Data", font = "Times 16 bold").grid(row= ROW, column=COL, sticky = W)
   COL+=1
   Button(fit_window, text="Plot Exp Data", command = plot_exp_data, font = "Times 16 bold").grid(row=ROW, column=COL, pady=4)
   COL-=1
   ROW += 1
   enter_num('fit_file', "Experimental Data Filename", fit_window,ROW, COL)
   ROW += 1
   enter_num('center', "Center of Beamstop (x y)", fit_window,ROW, COL)
   ROW += 1
   enter_num('border', "Additional Cropping", fit_window,ROW, COL)
   ROW += 1
   enter_num('mask_threshold', "Mask Threshhold", fit_window,ROW, COL)
   ROW += 1
   Label(fit_window, text="Fit Parameters", font = "Times 16 bold").grid(row= ROW, column=COL, columnspan =2, sticky = W)
   ROW += 1
   enter_num('background', "Background Noise", fit_window,ROW, COL)

   ROW+=1      # These are checkboxes which, if unchecked, will hold fixed fit parameters.
   tick("fit_radius_1", g.var_names[0], fit_window,ROW,COL)
   COL+=1
   tick("fit_radius_2", g.var_names[1], fit_window,ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_rho_1", g.var_names[3], fit_window,ROW,COL)
   COL+=1
   tick("fit_rho_2", g.var_names[4], fit_window,ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_z_dim", g.var_names[2], fit_window,ROW,COL)
   COL+=1
   tick("fit_x_theta", "x rotation", fit_window,ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_y_theta", "y rotation", fit_window,ROW,COL)
   COL+=1
   tick("fit_z_theta", "z rotation", fit_window,ROW,COL)
   COL-=1
   ROW+=1
   tick("fit_background", "background", fit_window,ROW,COL)
   COL+=1
   tick("fit_num", g.var_names[5], fit_window,ROW,COL)
   COL-=1
   ROW += 1
   tick("fit_length_2", g.var_names[6], fit_window,ROW,COL)
   ROW += 1
   for i in range(len(g.var_names)):
       if g.var_names[i] == "unused":
          g.dictionary_in['fit_'+g.var_list[i]+'2'].config(state=DISABLED)
       else:
          g.dictionary_in['fit_'+g.var_list[i]+'2'].config(state=NORMAL)
          g.labels[g.var_list[i]].config(state=NORMAL)
   
   enter_num('max_iter', "Maximum Iterations (0=default)", fit_window,ROW, COL)
   ROW += 1
   enter_num('update_freq', "Update Interval", fit_window,ROW, COL)   #TODO: debug
   ROW += 1
   enter_num('grid_compression', "Grid Compression (2, 5, or 10)", fit_window,ROW, COL)
   ROW += 1
   tick('plot_fit_tick',"Plot Fit Results", fit_window,ROW, COL)
   COL += 1
   tick('plot_residuals_tick',"Plot Fit Residuals", fit_window,ROW, COL)
   COL-=1
   ROW += 1
   Button(fit_window, text="Calc Residuals", command = plot_residuals, font = "Times 16 bold").grid(row=ROW, column=COL, pady=4)
   COL+= 1
   Button(fit_window, text="Fit Exp Data", command = perform_fit, font = "Times 16 bold").grid(row=ROW, column=COL, pady=4)
   COL -= 1


