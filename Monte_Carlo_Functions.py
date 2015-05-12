import random, os, sys, pylab, time
import numpy as np

global accelerated
global opencl_enabled,opencl_instance


######################          Finding the Point used in the Calculation         ##################

def Points_For_Calculation(seed=0):
    global dictionary_SI,quiet

    if seed:
       np.random.seed([seed])
    
    x_dim = dictionary_SI['x_dim']
    y_dim = dictionary_SI['y_dim']
    z_dim = dictionary_SI['z_dim']
    x_theta = dictionary_SI['x_theta']
    y_theta = dictionary_SI['y_theta']
    z_theta = dictionary_SI['z_theta']
    ave_dist = dictionary_SI['ave_dist']
    z_scale = dictionary_SI['z_scale']
    #I make a grid, then find a random number from a normal distribution with radius ave_dist. this gets added to the grid coordinates to randomise this.
    RandomPoints = np.asarray([((np.random.normal()*dictionary_SI['travel']+x_coord)%dictionary_SI['x_dim'] - dictionary_SI['x_dim']/2, (np.random.normal()*dictionary_SI['travel']+y_coord)%dictionary_SI['y_dim'] - dictionary_SI['y_dim']/2, (np.random.normal()*dictionary_SI['travel']*z_scale+z_coord)%dictionary_SI['z_dim'] - dictionary_SI['z_dim']/2)
                    for z_coord in np.arange(-z_dim/2, z_dim/2, ave_dist*z_scale) for y_coord in np.arange(-y_dim/2, y_dim/2, ave_dist) for x_coord in np.arange(-x_dim/2, x_dim/2, ave_dist)])


    points_inside = np.asarray([np.append(coords, [density(coords)])
                                  for coords in RandomPoints if abs(density(coords))>0.00001])
    
    #To use less RAM, i am clearing this variable now.
    RandomPoints = None

    if not quiet:
        print("{0} points will be used for the calculation.".format(len(points_inside)))
    sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\n"+str(len(points_inside)) + " points were used for the calculation.")
    sim_info.close()
    
    #Rotation Matricies. The last column and row of each matrix is so that you can keep the density the same when multiplying the matrix by each point.
    rotz = np.array([[np.cos(z_theta),-np.sin(z_theta),0,0],[np.sin(z_theta),np.cos(z_theta),0,0],[0,0,1,0],[0,0,0,1]])
    roty = np.array([[np.cos(y_theta),0,np.sin(y_theta),0],[0,1,0,0],[-np.sin(y_theta),0,np.cos(y_theta),0],[0,0,0,1]])
    rotx = np.array([[1,0,0,0],[0,np.cos(x_theta),-np.sin(x_theta),0],[0,np.sin(x_theta),np.cos(x_theta),0],[0,0,0,1]])

    #Multiplying the matrix by each column. The transpose is to make the multiplication work properely.
    #I am multiplying the rotation matricies together, then multiplying it by the coordinates for each point
    try:
       return np.asarray(points_inside.dot(np.transpose(rotz.dot(roty).dot(rotx))))
    except ValueError:
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
symmetric = dictionary_SI['symmetric']
Qz = dictionary_SI['Qz']

#for asymmetric objects, no small angle approximation
if opencl_enabled:
   qsize=dictionary_SI['qsize']
   ehc=dictionary_SI['ehc']
   pixels=dictionary_SI['pixels']
   return opencl_instance.sumint(qsize,ehc,pixels,Points,symmetric,Qz)
elif symmetric == 0 and Qz == 0:
    #print "No symmetry; no small angle approximation."
    def Detector_Intensity(Points,mask=[]):
        global dictionary_SI
        QSize = dictionary_SI['QSize']
        pixels = dictionary_SI['pixels']
        EHC = dictionary_SI['EHC']
        if accelerated:
            if not len(mask):
                mask = np.ones((pixels,pixels))
            #return fastmath.fastmath.sumintensity00(QSize,EHC,mask,Points)
            return fastmath.fastmath.sumintensity00(QSize,EHC,mask,Points.T)
        else:
            #print "FYI: Could not accelerate using f2py."
            Intensity = np.array([[np.sum(np.cos(np.sum(
                [row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*pixels)**2 + (col-0.5*pixels)**2)**0.5)*QSize/pixels/2/EHC)**2]
                *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
                          +np.sum(np.sin(np.sum(
                              [row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*pixels)**2 + (col-0.5*pixels)**2)**0.5)*QSize/pixels/2/EHC)**2]
                              *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
                          for col in range(int(pixels))] for row in range(int(pixels))])
            return Intensity/np.sum(Intensity)

#for asymmetric objects, small angle approximation
elif symmetric == 0 and Qz == 1:
    #print "No symmetry; small angle approximation."
    def Detector_Intensity(Points,mask=[]):
        global dictionary_SI
        QSize = dictionary_SI['QSize']
        pixels = dictionary_SI['pixels']
        EHC = dictionary_SI['EHC']
        if accelerated:
            if not len(mask):
                mask = np.ones((pixels,pixels))
            #return fastmath.fastmath.sumintensity00(QSize,EHC,mask,Points)    #Not a typo; sumint01 is (slightly) slower and thus pointless.
            return fastmath.fastmath.sumintensity00(QSize,EHC,mask,Points.T)  #Not a typo; sumint 01 is (slightly) slower and thus pointless.
        else:
            #print "FYI: Could not accelerate using f2py."
            Intensity = np.array([[np.sum(np.cos(np.sum(
                [row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize, 0]
                *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
                          +np.sum(np.sin(np.sum(
                              [row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*pixels)**2 + (col-0.5*pixels)**2)**0.5)*QSize/pixels/2/EHC)**2]
                              *Points[:,0:3], axis = 1))*np.transpose(Points[:,3:4]))**2
                          for col in range(int(pixels))] for row in range(int(pixels))])
            return Intensity/np.sum(Intensity)

#for symmetric objects, no small angle approximation
elif symmetric == 1 and Qz == 0:
    #print "Symmetry; no small angle approximation."
    def Detector_Intensity(Points,mask=[]):
        global dictionary_SI
        QSize = dictionary_SI['QSize']
        pixels = dictionary_SI['pixels']
        EHC = dictionary_SI['EHC']
        if accelerated:
            if not len(mask):
                mask = np.ones((pixels,pixels))
            #return fastmath.fastmath.sumintensity10(QSize,EHC,mask,Points)
            return fastmath.fastmath.sumintensity10(QSize,EHC,mask,Points.T)
        else:
            #print "FYI: Could not accelerate using f2py."
            Intensity = np.array([[np.sum(np.cos(np.sum(
                [row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize, 2*EHC*np.sin((((row-0.5*pixels)**2 + (col-0.5*pixels)**2)**0.5)*QSize/pixels/2/EHC)**2]
                *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
                          for col in range(int(pixels))] for row in range(int(pixels))])
            return Intensity/np.sum(Intensity)

#for symmetric objects, small angle approximation
elif symmetric == 1 and Qz == 1:
    #print "Symmetry; small angle approximation."
    def Detector_Intensity(Points,mask=[]):
        global dictionary_SI
        QSize = dictionary_SI['QSize']
        pixels = dictionary_SI['pixels']
        EHC = dictionary_SI['EHC']
        if accelerated:
            if not len(mask):
                mask = np.ones((pixels,pixels))
            #return fastmath.fastmath.sumintensity11(QSize,mask,Points)
            return fastmath.fastmath.sumintensity11(QSize,mask,Points.T)
        else:
            #print "FYI: Could not accelerate using f2py."
            Intensity = np.array([[np.sum(np.cos(np.sum(
                [row*QSize/pixels-0.5*QSize, col*QSize/pixels-0.5*QSize, 0.]
                *Points[:,0:3],axis =1))*np.transpose(Points[:,3:4]))**2
                          for col in range(int(pixels))] for row in range(int(pixels))])
            return Intensity/np.sum(Intensity)




###########          Average Intensity         #############

def Average_Intensity():
    global dictionary_SI,debug
    num_plots = dictionary_SI['num_plots']
    print "START TIME: "+time.strftime("%X")
    sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\nStart Time: "+time.strftime("%X"))
    sim_info.close()
    for plot_number in range(int(num_plots)):

        print "Average Plot " + str(plot_number+1) + " out of " + str(int(num_plots))
        
        sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
        sim_info.write("\nAverage Plot, plot " + str(plot_number+1) + " out of " + str(int(num_plots)) )
        sim_info.close()

        try:
            dictionary_SI['TEMP_VAR'] #This is here so it will only make the estimated time once.
            ##Intensity = Detector_Intensity(Points_For_Calculation())  #Commented and separated so I can time these separately.
            Points = Points_For_Calculation()
            Intensity = Detector_Intensity(Points)
            if debug:
               print("FINISHED CALCULATION {0}: {1}".format(plot_number+1,time.strftime("%X")))
        except KeyError:
            Points = Points_For_Calculation()
            try:
                dictionary_SI['current_value']
                est_time = 0# 10**-7*len(Points)*dictionary_SI['pixels']**2*dictionary_SI['num_plots']*dictionary_SI['s_step']
            except KeyError:
                est_time = 0# 10**-7*len(Points)*dictionary_SI['pixels']**2*dictionary_SI['num_plots']
            mins, secs = divmod(est_time, 60)
            hours, mins = divmod(est_time, 60)
            #print "Estimated time to finish all calculations: " + str(int(hours)) + " hours, " + str(int(mins)) + " minutes and " + str(int(secs)) + " seconds."
            dictionary_SI['TEMP_VAR'] = 0
            Intensity = Detector_Intensity(Points)
            if debug:
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
    sim_info = open(dictionary_SI['path_to_subfolder']+"simulation_infomation.txt","a")
    sim_info.write("\nEnd Time: "+time.strftime("%X"))
    sim_info.close()
    return cumulative/np.sum(cumulative)




###############################         Radial Intensity          #######################################
def radial(Intensity):
    global dictionary_SI
    QSize = dictionary_SI['QSize']
    pixels = dictionary_SI['pixels']
    EHC = dictionary_SI['EHC']
    num_plot_points = dictionary_SI['num_plot_points']
    delta = dictionary_SI['delta']
    return np.array([[0.5*QSize*temp/num_plot_points, np.mean([Intensity[x,y] for x in range(int(pixels)) for y in range(int(pixels))
                        if temp-0.5*delta <=  ( (x - 0.5*pixels)**2 + (y - 0.5*pixels)**2 )**0.5  <=temp+0.5*delta])]
                     for temp in range(int(num_plot_points))])[1:,]



#########################             Plot Angle at a fixed radius           ###########################
def plotting_circle(Intensity):
    global dictionary_SI
    pixels = dictionary_SI['pixels']
    xrow, yrow = np.shape(Intensity)
    
    circ_delta = float(dictionary_SI['circ_delta'])
    theta_delta = float(dictionary_SI['theta_delta'])
    pixel_radius = float(dictionary_SI['pixel_radius'])
    

    rad_theta = np.asarray([[  [ ((x-0.5*pixels)**2 + (y-0.5*pixels)**2)**0.5 , np.angle((x-0.5*pixels)+1j*(y-0.5*pixels))+3.1416  ] for y in range(int(yrow))] for x in range(int(xrow))])
    return np.array([[theta,
                      np.mean( [Intensity[x,y] for x in range(int(xrow))  for y in range(int(yrow))
                        if pixel_radius-circ_delta/2.         <=rad_theta[x,y,0]<=    pixel_radius+circ_delta/2.
                               and theta-theta_delta/2.       <=rad_theta[x,y,1]<=    theta+theta_delta/2. ])]
                                     for theta in np.arange(theta_delta/2. , 2.*3.1416-theta_delta/2. , theta_delta )])
    
        
                      
    

#################################           Saving Data to a CSV File           #####################
#Name must be entered as a string
def save(data, name):
    global dictionary_SI
    np.savetxt(dictionary_SI['path_to_subfolder']+name+".csv", data, delimiter=",")

