import random, os, sys, pylab, time
import numpy as np

import global_vars as g

######################          Finding the Point used in the Calculation         ##################

def Points_For_Calculation(seed=0):
    if seed:
       np.random.seed([seed])
    
    x_dim,y_dim,z_dim = g.dictionary_SI['x_dim'],g.dictionary_SI['y_dim'],g.dictionary_SI['z_dim']
    x_theta,y_theta,z_theta = g.dictionary_SI['x_theta'],g.dictionary_SI['y_theta'],g.dictionary_SI['z_theta']
    ave_dist = g.dictionary_SI['ave_dist']
    z_scale = g.dictionary_SI['z_scale']
    #I make a grid, then find a random number from a normal distribution with radius ave_dist. this gets added to the grid coordinates to randomise this.
    RandomPoints = np.asarray([((np.random.normal()*g.dictionary_SI['travel']+x_coord)%x_dim - x_dim/2, (np.random.normal()*g.dictionary_SI['travel']+y_coord)%y_dim - y_dim/2, (np.random.normal()*g.dictionary_SI['travel']*z_scale+z_coord)%z_dim - z_dim/2)
                    for z_coord in np.arange(-z_dim/2, z_dim/2, ave_dist*z_scale) for y_coord in np.arange(-y_dim/2, y_dim/2, ave_dist) for x_coord in np.arange(-x_dim/2, x_dim/2, ave_dist)])

    g.dprint("{0}: Generated {1} random points.".format(time.strftime("%X"),RandomPoints.shape[0]))

    #Fortran implementation about 8x faster than new python implementation.
    if g.accelerate_points and g.f2py_enabled and RandomPoints.shape[0] > 100000 and g.dictionary_SI['shape'] in (1,2,3,4,5,6,7,11,13,14,15,16):
     try:
        g.dprint('{0}: Using Fortran to calculate densities.'.format(time.strftime("%X")))
        densities = np.float32(np.append(RandomPoints,np.zeros([RandomPoints.shape[0],1]),1)).T
        if g.dictionary_SI['shape'] == 1:
            fastmath.density.d1sphere(g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1'],densities)
        elif g.dictionary_SI['shape'] == 2:
            fastmath.density.d2cylinder(g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1'],densities)
        elif g.dictionary_SI['shape'] == 3:
            fastmath.density.d3coreshell(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['rho_2'],densities)
        elif g.dictionary_SI['shape'] == 4:
            fastmath.density.d4gaussian(g.dictionary_SI['radius_2'],densities)
        elif g.dictionary_SI['shape'] == 5:
            fastmath.density.d5choppedcone(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['z_dim'],densities)
        elif g.dictionary_SI['shape'] == 6:
            fastmath.density.d6hexprism(g.dictionary_SI['radius_1'],g.dictionary_SI['rho_1'],densities)
        elif g.dictionary_SI['shape'] == 7:
            fastmath.density.d7rectprism(g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],densities)
        elif g.dictionary_SI['shape'] == 11:
            fastmath.density.d11doubleslit(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],densities)
        elif g.dictionary_SI['shape'] == 13:
            fastmath.density.d13sine(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['rho_2'],g.dictionary_SI['z_dim'],densities)
        elif g.dictionary_SI['shape'] == 14:
            fastmath.density.d14doublecone(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['z_dim'],densities)
        elif g.dictionary_SI['shape'] == 15:
            fastmath.density.d15elipticalcylinder(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],densities)
        elif g.dictionary_SI['shape'] == 16:
            fastmath.density.d16asymmhexpyr(g.dictionary_SI['radius_1'],g.dictionary_SI['radius_2'],g.dictionary_SI['rho_1'],g.dictionary_SI['z_dim'],densities)
        densities = densities.T
        outside = [i for i in range(densities.shape[0]) if not densities[i,3]]
        points_inside = np.delete(densities,outside,axis=0)
     except AttributeError:    #In case fortran binary is too old.
        print("Could not speed up with fortran.  Recompile.")
        points = np.c_[RandomPoints,density(RandomPoints)]
        outside = [i for i in range(points.shape[0]) if not points[i,3]]
        points_inside = np.delete(points,outside,axis=0)
    elif g.accelerate_points and g.opencl_enabled and RandomPoints.shape[0] > 100000 and g.dictionary_SI['shape'] in (1,2,3,4,5,6,7,11,13,14):
        g.dprint('{0}: Using OpenCL for density calculation.'.format(time.strftime("%X")))
        densities = g.opencl_density.density(RandomPoints)
        points = np.c_[RandomPoints,densities]
        outside = [i for i in range(points.shape[0]) if not points[i,3]]
        points_inside = np.delete(points,outside,axis=0)
    else:
        points = np.c_[RandomPoints,density(RandomPoints)]
        outside = [i for i in range(points.shape[0]) if not points[i,3]]
        points_inside = np.delete(points,outside,axis=0)
        #points_inside = np.asarray([np.append(coords, [density(coords)]) for coords in RandomPoints if abs(density(coords))>0.00001])  #30% slower implementation
    
    RandomPoints = None         #To use less RAM, i am clearing this variable now.

    if not g.quiet:
        print("{0}: {1} points will be used for the calculation.".format(time.strftime("%X"),len(points_inside)))
    sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\n"+str(len(points_inside)) + " points were used for the calculation.")
    sim_info.close()
    
    #Rotation Matricies. The last column and row of each matrix is so that you can keep the density the same when multiplying the matrix by each point.
    rotz = np.array([[np.cos(z_theta),-np.sin(z_theta),0,0],[np.sin(z_theta),np.cos(z_theta),0,0],[0,0,1,0],[0,0,0,1]])
    roty = np.array([[np.cos(y_theta),0,np.sin(y_theta),0],[0,1,0,0],[-np.sin(y_theta),0,np.cos(y_theta),0],[0,0,0,1]])
    rotx = np.array([[1,0,0,0],[0,np.cos(x_theta),-np.sin(x_theta),0],[0,np.sin(x_theta),np.cos(x_theta),0],[0,0,0,1]])

    #Multiplying the matrix by each column. The transpose is to make the multiplication work properly.
    #I am multiplying the rotation matricies together, then multiplying it by the coordinates for each point
    try:
       return np.asarray(points_inside.dot(np.transpose(rotz.dot(roty).dot(rotx))))
       #return np.asarray(points_inside[points_inside[:,2].argsort()].dot(np.transpose(rotz.dot(roty).dot(rotx))))   #first orders list by z and then rotation matrices....backwards?
    except ValueError:
       print points.shape
       print rotx.shape, roty.shape, rotz.shape
       print rotz.dot(roty).dot(rotx).shape
       print points_inside.shape
       print np.transpose(rotz.dot(roty).dot(rotx)).shape
       print points_inside.dot(np.transpose(rotz.dot(roty).dot(rotx))).shape
       return np.asarray(points_inside.dot(np.transpose(rotz.dot(roty).dot(rotx))))

    #return np.asarray(points_inside.dot(np.transpose(rotz.dot(roty).dot(rotx))))





##################         Finding the Detector Intensity        #########################
#NOTES:
#Intensity = [[np.sum(np.cos(np.sum(<- (1,2,3)*(1,2,3)=(1,4,9) This sum adds the three components together to make this a dot product.
#For the dot product, Qz =  2*EHC*sin( sqrt(x**2 + y**2) *QSize/pixels/2/EHC )**2
#symmetric = g.dictionary_SI['symmetric']
#Qz = g.dictionary_SI['Qz']

def Calculate_Intensity(Points,mask=[],coherence_dup = 1, coherence_taper = 0):
   '''Runs Detector_Intensity, but can also handle objects longer than the coherence length.'''
   x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
   Points = Points[Points[:,2].argsort()]    #orders by z
   z_list = Points[:,2]
   #print(points_inside[points_inside[:,2].argsort()][::100,2])
   length = z_list[-1]-z_list[0]
   #coherence_length = 5e-7
   if not g.dictionary['d_lambda']:
      coherence_length = 1
   else:
      coherence_length = 2*np.pi/(g.dictionary_SI['EHC'] * g.dictionary['d_lambda'])
   #coherence_dup = 1 #number of bunches per coherence_length
   num_bunches = int(round(coherence_dup*length/coherence_length))
   if num_bunches > coherence_dup:    #if length > coherence_length:
      print("Object length ({0}) exceeds coherence length ({1})...".format(length,coherence_length))
      print("Will divide into {0} sections.".format(num_bunches))
      dividing_points = np.searchsorted(z_list,-length/2+(np.arange(num_bunches+1))*coherence_length/coherence_dup)
      dividing_points[0] = 0    #Not sure why this isn't already 0.
      dividing_points[-1] = len(z_list)  #Should this be the last element or should I append?
      if g.debug:
         print("Dividing points are:")
         print(dividing_points)
         #print([z_list[i] for i in dividing_points])
         #print([z_list[dividing_points[i+1]]-z_list[dividing_points[i] for i in range(len(dividing_points)-1)])
      intensity = np.zeros((y_pixels,x_pixels),order='F')   #todo: might not need Fortran order
      #TODO: DO I NEED TO OFFSET Z IN SOME WAY??
      if coherence_taper:
          for i in range(1,coherence_dup):
             if g.debug and g.verbose:
                print("Starting section {0} of {1}".format(i,len(dividing_points)+coherence_dup-1))
                print("{0} to {1}".format(dividing_points[0],dividing_points[i]))
                print("{0} to {1}".format(z_list[dividing_points[0]],z_list[dividing_points[i]]))
             intensity += Detector_Intensity(Points[dividing_points[0]:dividing_points[i],:],mask)
      for i in range(len(dividing_points)-coherence_dup):     #TODO: also do more at ends??
         if g.debug and g.verbose:
            print("Starting main section {0} of {1}".format(i+coherence_dup,len(dividing_points)+coherence_dup-1))
            #print("Starting section {0} of {1}".format(i+1,len(dividing_points)-coherence_dup))
            print("{0} to {1}".format(dividing_points[i],dividing_points[i+coherence_dup]-1))
            print("{0} to {1}".format(z_list[dividing_points[i]],z_list[dividing_points[i+coherence_dup]-1]))
         intensity += Detector_Intensity(Points[dividing_points[i]:dividing_points[i+coherence_dup],:],mask)
      if coherence_taper:
          for i in range(coherence_dup-1):
             if g.debug and g.verbose:
                print("Starting section {0} of {1}".format(i+1+len(dividing_points),len(dividing_points)+coherence_dup-1))
                print("{0} to {1}".format(dividing_points[len(dividing_points)-coherence_dup+i],dividing_points[-1]))
                print("{0} to {1}".format(z_list[dividing_points[len(dividing_points)-coherence_dup+i]],z_list[-1]))
             intensity += Detector_Intensity(Points[dividing_points[len(dividing_points)-coherence_dup+i]:dividing_points[-1],:],mask)
      return intensity
   else:
      if coherence_length < 1:
         g.vprint("Object length ({0}) is less than coherence length ({1})...".format(length,coherence_length))
      return Detector_Intensity(Points,mask)

if g.opencl_enabled:
   def Detector_Intensity(Points,mask=[]):
      symmetric = g.dictionary_SI['symmetric']
      qsize=g.dictionary_SI['QSize']
      ehc=g.dictionary_SI['EHC']
      x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
      #g.dprint("Shape {0}".format(g.dictionary['shape']))
      if symmetric:
         g.dprint("Using symmetry.")
      if not len(mask) or g.dictionary['grid_compression'] < 2:
         return g.opencl_sumint.sumint(qsize,ehc,x_pixels,y_pixels,Points,symmetric)
      else:
         if g.dictionary['grid_compression'] < 5:
            print('Grid compression of < 5 does not produce significant speedup when using OpenCL.  Using 0/1 or 5 or 10 is recommended.')
         return g.opencl_sumint.sumint_mask(qsize,ehc,mask,Points,symmetric)

elif g.f2py_enabled:
   def Detector_Intensity(Points,mask=[]):
      symmetric = g.dictionary_SI['symmetric']
      QSize = g.dictionary_SI['QSize']
      x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
      EHC = g.dictionary_SI['EHC']
      #g.dprint("Shape {0}".format(g.dictionary['shape']))
      if not len(mask):
         mask = np.ones((y_pixels,x_pixels))
      if symmetric:
         g.dprint("Using symmetry.")
         return fastmath.sumint.symmetric(QSize,mask,Points.T)
      else:
         return fastmath.sumint.asymmetric(QSize,EHC,mask,Points.T)

else:   #python only
   def Detector_Intensity(Points,mask=[]):
      symmetric = g.dictionary_SI['symmetric']
      QSize = g.dictionary_SI['QSize']
      x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
      EHC = g.dictionary_SI['EHC']
      #g.dprint("Shape {0}".format(g.dictionary['shape']))
      if symmetric:
         g.dprint("Using symmetry.")
         Intensity = np.array([[np.sum(np.cos(np.sum(
               [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
               *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
                         for col in range(int(x_pixels))] for row in range(int(y_pixels))])
         return Intensity/np.sum(Intensity)
      else:
         Intensity = np.array([[np.sum(np.cos(np.sum(
            [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
              *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
                        +np.sum(np.sin(np.sum(
                            [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
                            *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
                        for col in range(int(x_pixels))] for row in range(int(y_pixels))])
         return Intensity/np.sum(Intensity)



##for asymmetric objects, no small angle approximation
#elif symmetric == 0 and Qz == 0:
#    #print "No symmetry; no small angle approximation."
#    def Detector_Intensity(Points,mask=[]):
#        QSize = g.dictionary_SI['QSize']
#        x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
#        EHC = g.dictionary_SI['EHC']
#        if g.f2py_enabled:
#            if not len(mask):
#                mask = np.ones((y_pixels,x_pixels))
#            return fastmath.sumint.sumintensity00(QSize,EHC,mask,Points.T)
#        else:
#            Intensity = np.array([[np.sum(np.cos(np.sum(
#                [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
#                *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
#                          +np.sum(np.sin(np.sum(
#                              [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
#                              *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
#                          for col in range(int(x_pixels))] for row in range(int(y_pixels))])
#            return Intensity/np.sum(Intensity)
#
##for asymmetric objects, small angle approximation
#elif symmetric == 0 and Qz == 1:
#    #print "No symmetry; small angle approximation."
#    def Detector_Intensity(Points,mask=[]):
#        QSize = g.dictionary_SI['QSize']
#        x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
#        EHC = g.dictionary_SI['EHC']
#        if g.f2py_enabled:
#            if not len(mask):
#                mask = np.ones((y_pixels,x_pixels))
#            return fastmath.sumint.sumintensity00(QSize,EHC,mask,Points.T)  #Not a typo; sumint 01 is (slightly) slower and thus pointless.
#        else:
#            Intensity = np.array([[np.sum(np.cos(np.sum(
#                [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 0]
#                *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
#                          +np.sum(np.sin(np.sum(
#                              [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
#                              *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
#                          for col in range(int(x_pixels))] for row in range(int(y_pixels))])
#            return Intensity/np.sum(Intensity)
#
#
##for symmetric objects, no small angle approximation
#elif symmetric == 1 and Qz == 0:
#    #print "Symmetry; no small angle approximation."
#    def Detector_Intensity(Points,mask=[]):
#        QSize = g.dictionary_SI['QSize']
#        x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
#        EHC = g.dictionary_SI['EHC']
#        if g.f2py_enabled:
#            if not len(mask):
#                mask = np.ones((y_pixels,x_pixels))
#            return fastmath.sumint.sumintensity10(QSize,EHC,mask,Points.T)
#        else:
#            Intensity = np.array([[np.sum(np.cos(np.sum(
#                [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/x_pixels/2/EHC)**2]
#                *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
#                          for col in range(int(x_pixels))] for row in range(int(y_pixels))])
#            return Intensity/np.sum(Intensity)
#
##for symmetric objects, small angle approximation
#elif symmetric == 1 and Qz == 1:
#    #print "Symmetry; small angle approximation."
#    def Detector_Intensity(Points,mask=[]):
#        QSize = g.dictionary_SI['QSize']
#        x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
#        EHC = g.dictionary_SI['EHC']
#        if g.f2py_enabled:
#            if not len(mask):
#                mask = np.ones((y_pixels,x_pixels))
#            return fastmath.sumint.sumintensity11(QSize,mask,Points.T)
#        else:
#            Intensity = np.array([[np.sum(np.cos(np.sum(
#                [row*QSize/y_pixels-0.5*QSize, col*QSize/x_pixels-0.5*QSize, 0.]
#                *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
#                          for col in range(int(x_pixels))] for row in range(int(y_pixels))])
#            return Intensity/np.sum(Intensity)




###########          Average Intensity         #############

def Average_Intensity(mask=[]):
    num_plots = g.dictionary_SI['num_plots']
    print "START TIME: "+time.strftime("%X")
    sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\nStart Time: "+time.strftime("%X"))
    sim_info.close()
    for plot_number in range(int(num_plots)):

        print "Average Plot " + str(plot_number+1) + " out of " + str(int(num_plots))
        
        sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
        sim_info.write("\nAverage Plot, plot " + str(plot_number+1) + " out of " + str(int(num_plots)) )
        sim_info.close()

        try:
            g.dictionary_SI['TEMP_VAR'] #This is here so it will only make the estimated time once.
            ##Intensity = Calculate_Intensity(Points_For_Calculation())  #Commented and separated so I can time these separately.
            Points = Points_For_Calculation()
            Intensity = Calculate_Intensity(Points,mask)
            g.vprint("FINISHED CALCULATION {0}: {1}".format(plot_number+1,time.strftime("%X")))
        except KeyError:
            Points = Points_For_Calculation()
            try:
                g.dictionary_SI['current_value']
                est_time = 0# 10**-7*len(Points)*g.dictionary_SI['pixels']**2*g.dictionary_SI['num_plots']*g.dictionary_SI['s_step']
            except KeyError:
                est_time = 0# 10**-7*len(Points)*g.dictionary_SI['pixels']**2*g.dictionary_SI['num_plots']
            mins, secs = divmod(est_time, 60)
            hours, mins = divmod(est_time, 60)
            #print "Estimated time to finish all calculations: " + str(int(hours)) + " hours, " + str(int(mins)) + " minutes and " + str(int(secs)) + " seconds."
            g.dictionary_SI['TEMP_VAR'] = 0
            Intensity = Calculate_Intensity(Points,mask)
            if g.verbose > 0:
               print("FINISHED CALCULATION {0}: {1}".format(plot_number+1,time.strftime("%X")))
            else:
               print "FINISHED FIRST CALCULATION: "+time.strftime("%X")

        
        try:
            cumulative += np.asarray(Intensity)
            #this finds the cumulative intensity
            #the "try" needs to be here as I have not defined cumulative from the start.
        except NameError:
            cumulative = np.asarray(Intensity)
    #this finds the average of all plots
    print "END TIME: "+time.strftime("%X")
    sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\nEnd Time: "+time.strftime("%X"))
    sim_info.close()
    return cumulative/np.sum(cumulative)




###############################         Radial Intensity          #######################################
def radial(Intensity):
    QSize = g.dictionary_SI['QSize']
    #x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
    x_pixels,y_pixels=Intensity.shape
    EHC = g.dictionary_SI['EHC']
    num_plot_points = g.dictionary_SI['num_plot_points']
    delta = g.dictionary_SI['delta']
    return np.array([[0.5*QSize*temp/num_plot_points, np.mean([Intensity[x,y] for x in range(x_pixels) for y in range(y_pixels) if temp-0.5*delta <=  np.sqrt( (x - 0.5*x_pixels)**2 + (y - 0.5*y_pixels)**2 )  <=temp+0.5*delta])] for temp in range(int(num_plot_points))])[1:,]



#########################             Plot Angle at a fixed radius           ###########################
def plotting_circle(Intensity):
    x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
    xrow, yrow = np.shape(Intensity)
    
    circ_delta = float(g.dictionary_SI['circ_delta'])
    theta_delta = float(g.dictionary_SI['theta_delta'])
    pixel_radius = float(g.dictionary_SI['pixel_radius'])
    

    rad_theta = np.asarray([[  [ ((x-0.5*x_pixels)**2 + (y-0.5*y_pixels)**2)**0.5 , np.angle((x-0.5*x_pixels)+1j*(y-0.5*y_pixels))+3.1416  ] for y in range(int(yrow))] for x in range(int(xrow))])
    return np.array([[theta,
                      np.mean( [Intensity[x,y] for x in range(int(xrow))  for y in range(int(yrow))
                        if pixel_radius-circ_delta/2.         <=rad_theta[x,y,0]<=    pixel_radius+circ_delta/2.
                               and theta-theta_delta/2.       <=rad_theta[x,y,1]<=    theta+theta_delta/2. ])]
                                     for theta in np.arange(theta_delta/2. , 2.*3.1416-theta_delta/2. , theta_delta )])
    
        
                      
    

#################################           Saving Data to a CSV File           #####################
#Name must be entered as a string
def save(data, name):
    np.savetxt(g.dictionary_SI['path_to_subfolder']+name+".csv", data, delimiter=",")

