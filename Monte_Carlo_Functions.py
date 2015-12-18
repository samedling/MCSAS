import random, os, sys, pylab, time
import numpy as np
from scipy.interpolate import interp2d,Rbf,griddata

import global_vars as g

######################          Finding the Point used in the Calculation         ##################

def rot_points(points,reverse=True):
   '''Rotates, or counter-rotates, a list of points.'''
   x_theta=g.dictionary_SI['x_theta']
   y_theta=g.dictionary_SI['y_theta']
   z_theta=g.dictionary_SI['z_theta']
   rotx=np.array([[1,0,0,0],[0,np.cos(x_theta),-np.sin(x_theta),0],[0,np.sin(x_theta),np.cos(x_theta),0],[0,0,0,1]])
   roty=np.array([[np.cos(y_theta),0,np.sin(y_theta),0],[0,1,0,0],[-np.sin(y_theta),0,np.cos(y_theta),0],[0,0,0,1]])
   rotz=np.array([[np.cos(z_theta),-np.sin(z_theta),0,0],[np.sin(z_theta),np.cos(z_theta),0,0],[0,0,1,0],[0,0,0,1]])   
   if reverse:
      return points.dot(rotz.dot(roty).dot(rotx))
   else:
      return points.dot(np.transpose(rotz.dot(roty).dot(rotx)))



def Points_For_Calculation(seed=0,sort=0):
    if seed:
       np.random.seed([seed])
    
    density(np.asarray([[1,1,1],[2,2,2]])) #This is so that you can manually redefine x_dim and y_dim in the density function. e.g. for the gaussian model, you may want to make more points.
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
       if sort:
          return np.asarray(points_inside[points_inside[:,2].argsort()].dot(np.transpose(rotz.dot(roty).dot(rotx))))    #first orders points by z, then rotates
       else:
          return np.asarray(points_inside.dot(np.transpose(rotz.dot(roty).dot(rotx))))
    except ValueError:
       if len(points_inside) == 0:
          print("Error: shape has no points; check your parameters.")
       else:
          print("Unknown error: it might work if you just run it again without changing anything.  Sorry.")
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


def Accurate_Intensity(Points,mask=[],sort=0):
   '''SLOW but accurately accounts for coherence length.'''
   symmetric = g.dictionary_SI['symmetric']
   QSize = g.dictionary_SI['QSize']
   x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
   EHC = g.dictionary_SI['EHC']
   coherence_length = 2*np.pi/(g.dictionary_SI['EHC']*g.dictionary['d_lambda'])
   #print('Object length is {0}'.format(Points[-1,2]-Points[0,2])) #wrong; points not sorted
   #print('Coherence length is {0}'.format(coherence_length))
   if g.f2py_enabled:  #fortran
      if not len(mask):
         mask = np.ones((y_pixels,x_pixels))    #TODO: Is this right or backwards?
      g.dprint('Calling Fortran...')
      if sort:
         return fastmath.sumint.coherent_shorter(QSize,EHC,coherence_length,mask,Points.T)
      else:
         return fastmath.sumint.coherent_long(QSize,EHC,coherence_length,mask,Points.T)
   elif g.opencl_enabled:   #OpenCL
      g.dprint('Calling OpenCL...')
      if not len(mask):
         return g.opencl_sumint.sumint(QSize,EHC,x_pixels,y_pixels,Points,symmetric,coherence_length=coherence_length)
      else:
         return g.opencl_sumint.sumint_mask(QSize,EHC,mask,Points,symmetric,coherence_length=coherence_length)
   else:                      #python
      print('Using Python (will probably take forever, and not debugged!)...')
      for i in range(x_pixels):
         for j in range(y_pixels):
            if (mask[j,i] > 0):
               Q = [i*QSize/x_pixels-0.5*QSize, j*QSize/y_pixels-0.5*QSize,
               2*EHC*np.sin(np.sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*QSize/(y_pixels*2*EHC))**2 ]
               temp_intensity = 0
               for p1 in range(len(Points)):
                  for p2 in range(len(Points)):
                     r = Points[p1,:] - Points[p2,:]
                     if (np.sum(r**2) > coherence_length**2):
                        QdotR = np.dot(Q,r)
                        temp_intensity += Points[p1,3]*Points[p2,3]*np.cos(QdotR)
      return intensity

def Adjust_Intensity(points,mask=[],newmask_shape=(20,20),newmask_points=400,interp_method='rbf'):
    '''Runs Accurate_Intensity on a small subset of the points in order to correct the faster Calculate_Intensity.'''
    #Create New Mask
    x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
    if not len(mask):
        newmask = np.zeros((x_pixels,y_pixels)) #TODO: Is this right or backwards?
        modx = x_pixels/newmask_shape[0]
        mody = y_pixels/newmask_shape[1]
        for i in range(newmask_shape[0]):
            for j in range(newmask_shape[1]):
                newmask[int((i+0.5)*modx),int((j+0.5)*mody)] = 1
    else:
        newmask = mask.copy()
        mask_points = mask.sum()
        mod = mask_points/newpoints
        counter = int((mask_points%newpoints)/2)
        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                if mask[i,j]:
                    counter = (counter+1)%mod
                    if counter == 1:
                        newmask[i,j] = 0
    #Calculate Accurate Intensity at New Mask points, and Fast Intensity at all points
    g.dprint('Calculating actual intensity on {0} pixels...'.format(int(newmask.sum())))
    accurate = Accurate_Intensity(points,newmask)
    g.dprint('Calculating rough intensity everywhere...')
    fast = Calculate_Intensity(points,mask)
    scaleby = normalize(accurate)/normalize(fast)
    
    #Interpolation Options: scipy.interpolate: Rbf, interp2d, or griddata
    if interp_method[0] == 'r':
        #Rbf
        x = []
        y = []
        z = []
        for i in range(x_pixels):
            for j in range(y_pixels):
                if newmask[i,j]:
                    x.append(i)
                    y.append(j)
                    z.append(scaleby[i,j])
        xx,yy = np.meshgrid(np.arange(x_pixels),np.arange(y_pixels))
        rbf = Rbf(x,y,z,epsilon=2)
        scaleby = rbf(xx,yy)
    
    #interp2d
    elif interp_method[0] == 'i':
        x = []
        y = []
        z = []
        for i in range(x_pixels):
            for j in range(y_pixels):
                if newmask[i,j]:
                    x.append(i)
                    y.append(j)
                    z.append(scaleby[i,j])
        xx,yy = np.meshgrid(x,y)
        f = interp2d(x,y,z,kind='cubic')    # 'cubic', 'linear', or 'nearest'
        newx = np.arange(x_pixels)
        newy = np.arange(y_pixels)
        scaleby = f(newx,newy)
    
    #griddata - seems broken; see http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html ??
    elif interp_method[0] == 'g':
        xy = []
        z = []
        for i in range(x_pixels):
            for j in range(y_pixels):
                if newmask[i,j]:
                    xy.append([i,j])
                    z.append(scaleby[i,j])
        xy = np.asarray(xy)
        gridx,gridy=np.mgrid[0:x_pixels,0:y_pixels]
        scaleby = griddata(xy,z,(gridx,gridy),method='cubic')
    
    print('Done Interpolating.')
    
    return fast*scaleby

def Calculate_Intensity(Points,mask=[],coherence_dup = 1, coherence_taper = 0):
   '''Runs Detector_Intensity, but can also handle objects longer than the coherence length with reasonable speed and accuracy.'''
   x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
   #Points = Points[Points[:,2].argsort()]    #orders by z
   #z_list = Points[:,2]   #WRONG if things have been rotated.
   z_list = rot_points(Points,True)[:,2]  #un-rotates
   Points = Points[z_list.argsort()]      #orders by unrotated z
   z_list = z_list[z_list.argsort()]      #orders by unrotated z
   #print(points_inside[points_inside[:,2].argsort()][::100,2])
   length = z_list[-1]-z_list[0]
   if len(z_list) == 0:
      print("Error: shape has no points; check your parameters.")
      return
   elif len(z_list) < 100:
      print length
      print("Warning: shape has few points; check your parameters.")
   #coherence_length = 5e-7
   if not g.dictionary['d_lambda']:
      coherence_length = 1
   else:
      coherence_length = 2*np.pi/(g.dictionary_SI['EHC'] * g.dictionary['d_lambda'])
   #num_bunches = int(round(coherence_dup*length/coherence_length))
   num_bunches = int(np.floor(coherence_dup*length/coherence_length))   #Seems best for some reason.
   #num_bunches = int(np.ceil(coherence_dup*length/coherence_length))   #worse than round, at least sometimes
   if num_bunches > coherence_dup:    #if length > coherence_length:
      piece_length = g.dictionary_SI['z_dim']/num_bunches
      print("Object length ({0:4.4} nm) exceeds coherence length ({1:4.4} nm)...".format(length*10**9,coherence_length*10**9))
      print("Will divide into {0} sections of length {1}.".format(num_bunches,piece_length))
      dividing_points = np.searchsorted(z_list,-length/2+(np.arange(num_bunches+1))*piece_length/coherence_dup)
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
            print("Starting main section {0} of {1}".format(i+coherence_dup,len(dividing_points)+coherence_dup-2))
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
         g.vprint("Object length ({0:4.4} nm) is less than coherence length ({1:4.4} nm)...".format(length*10**9,coherence_length*10**9))
      return Detector_Intensity(Points,mask)

if g.opencl_enabled:
   def Detector_Intensity(Points,mask=[]):
      symmetric = g.dictionary_SI['symmetric']
      Qz = g.dictionary_SI['Qz']
      qsize=g.dictionary_SI['QSize']
      ehc=g.dictionary_SI['EHC']
      x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
      #g.dprint("Shape {0}".format(g.dictionary['shape']))
      if not len(mask) or g.dictionary['grid_compression'] < 2:
         return g.opencl_sumint.sumint(qsize,ehc,x_pixels,y_pixels,Points,symmetric)
      else:
         if g.dictionary['grid_compression'] < 5:
            print('Grid compression of < 5 does not produce significant speedup when using OpenCL.  Using 0/1 or 5 or 10 is recommended.')
         return g.opencl_sumint.sumint_mask(qsize,ehc,mask,Points,symmetric)

elif g.f2py_enabled:
   def Detector_Intensity(Points,mask=[]):
      symmetric = g.dictionary_SI['symmetric']
      Qz = g.dictionary_SI['Qz']
      QSize = g.dictionary_SI['QSize']
      x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
      EHC = g.dictionary_SI['EHC']
      #g.dprint("Shape {0}".format(g.dictionary['shape']))
      if not len(mask):
         mask = np.ones((y_pixels,x_pixels))
      if symmetric:
         if Qz:
            g.dprint("Using both radial symmetry and small angle approximation.")
            return fastmath.sumint.symmetric_small(QSize,EHC,mask,Points.T)
         else:
            g.dprint("Using radial symmetry but not small angle approximation.")
            return fastmath.sumint.symmetric(QSize,EHC,mask,Points.T)
      else:
         if Qz:
            g.dprint("Using small angle approximation but not radial symmetry.")
            return fastmath.sumint.asymmetric_small(QSize,EHC,mask,Points.T)
         else:
            g.dprint("Not using radial symmetry or small angle approximation.")
            return fastmath.sumint.asymmetric(QSize,EHC,mask,Points.T)

else:   #python only
   def Detector_Intensity(Points,mask=[]):
      symmetric = g.dictionary_SI['symmetric']
      Qz = g.dictionary_SI['Qz']
      QSize = g.dictionary_SI['QSize']
      x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
      max_pixels=max(x_pixels,y_pixels)
      EHC = g.dictionary_SI['EHC']
      #g.dprint("Shape {0}".format(g.dictionary['shape']))
      if symmetric:
         if Qz:
            g.dprint("Using both radial symmetry and small angle approximation.")
            Intensity = np.array([[np.sum(np.cos(np.sum(
               [(row-0.5*y_pixels)*QSize/max_pixels, (col-0.5*x_pixels)*QSize/max_pixels, 0.]
               *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
               for col in range(int(x_pixels))] for row in range(int(y_pixels))])
         else:
            g.dprint("Using radial symmetry but not small angle approximation.")
            Intensity = np.array([[np.sum(np.cos(np.sum(
               [(row-0.5*y_pixels)*QSize/max_pixels, (col-0.5*x_pixels)*QSize/max_pixels, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/max_pixels/2/EHC)**2]
               *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
               for col in range(int(x_pixels))] for row in range(int(y_pixels))])
      else:
         if Qz:
            g.dprint("Using small angle approximation but not radial symmetry.")
            Intensity = np.array([[np.sum(np.cos(np.sum(
               [(row-0.5*y_pixels)*QSize/max_pixels, (col-0.5*x_pixels)*QSize/max_pixels, 0.]
               *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
               +np.sum(np.sin(np.sum(
               [(row-0.5*y_pixels)*QSize/max_pixels, (col-0.5*x_pixels)*QSize/max_pixels, 0.]
               *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
               for col in range(int(x_pixels))] for row in range(int(y_pixels))])
         else:
            g.dprint("Not using radial symmetry or small angle approximation.")
            Intensity = np.array([[np.sum(np.cos(np.sum(
               [(row-0.5*y_pixels)*QSize/max_pixels, (col-0.5*x_pixels)*QSize/max_pixels, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/max_pixels/2/EHC)**2]
               *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
               +np.sum(np.sin(np.sum(
               [(row-0.5*y_pixels)*QSize/max_pixels, (col-0.5*x_pixels)*QSize/max_pixels, 2*EHC*np.sin((((row-0.5*y_pixels)**2 + (col-0.5*x_pixels)**2)**0.5)*QSize/max_pixels/2/EHC)**2]
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
            print Points
            Intensity = Calculate_Intensity(Points,mask)
            g.vprint("FINISHED CALCULATION {0}: {1}".format(plot_number+1,time.strftime("%X")))
        except KeyError:
            Points = Points_For_Calculation()
            print Points
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
               print("FINISHED FIRST CALCULATION: "+time.strftime("%X"))

        
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


######################        Interparticle Scattering      #############################
def inter_intensity(mask=[]):
    print "START TIME: "+time.strftime("%X")
    sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\nStart Time: "+time.strftime("%X"))
    sim_info.write("\nx centre,   ycentre")
    AllPoints = None
    for temp in range(int(g.dictionary_SI['numinter'])):

        xcentre=np.random.random()*g.dictionary_SI['xinter']*10**-9
	ycentre=np.random.random()*g.dictionary_SI['yinter']*10**-9
        sim_info.write("\n"+str(xcentre) + " ,    " + str(ycentre))
        print xcentre,ycentre
        TempPoints = np.array([x+[xcentre,ycentre,0,0] for x in Points_For_Calculation()])
        try:
	        overlapdensity = density(np.array( [x[0:3]-[xcentre,ycentre,0] for x in AllPoints]))
		points = np.c_[AllPoints, overlapdensity]
		inside = [i for i in range(points.shape[0]) if points [i,4] ]
		points_to_keep = np.delete(AllPoints,inside,axis=0)
		AllPoints = np.array(np.append(TempPoints,points_to_keep, axis=0))
        except TypeError:
		AllPoints = np.array(TempPoints)
    sim_info.write("\nTotal Points Used in Calculation:"+str(len(AllPoints)))
    sim_info.close()
    #Points_Plot(AllPoints,'points',1)
    Intensity = Calculate_Intensity(AllPoints,mask)

    print "END TIME: "+time.strftime("%X")
    sim_info = open(g.dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\nEnd Time: "+time.strftime("%X"))
    sim_info.close()
    return Intensity/np.sum(Intensity)




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

