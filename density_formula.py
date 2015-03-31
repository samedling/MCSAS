#THIS IS THE MOST IMPORTANT LINE!!!!!!!!!!!!!!!!
#When making points, x<x_dim/2, y<y_dim/2. When defining a function, by default, the x dimensions and y dimensions are equal to double radius_1, i.e. 2*radius_1 = x_dim = y_dim.
#Make sure this holds true for your function.
#you can always override x_dim and y_dim by defining them in your program (e.g. dictionary_SI['x_dim']=5*10**-9)
#This will NOT!! work for a sequence though.

#A reminder: the camera is a long way away in the z direction
#The function is run once for each point, which is entered as a numpy array, i.e. it is in the form: np.array([x,y,z])
#Make sure you always use dictionary_SI to get variables.



import pylab, sys, os, random
import numpy as np

if dictionary_SI == 0:
    print "Analytic Model Only"

elif dictionary_SI['shape'] == 1:
    print "Monte Carlo Sphere"
    def density(coords):
        dist = np.sqrt(np.sum(coords**2))
        if dist < dictionary_SI['radius_1']:
            b=dictionary['rho_1']
        else:
            b=0
        return b

elif dictionary_SI['shape'] == 2:
    print "Monte Carlo Cylinder"
    def density(coords):
        dist = np.sqrt(np.sum(coords[0:2]**2))
        if dist < dictionary_SI['radius_1']:
            b=dictionary_SI['rho_1']
        else:
            b=0
        return b

elif dictionary_SI['shape'] == 3:
    print "Monte Carlo core shell cylinder"
    def density(coords):
        dist = np.sqrt(np.sum(coords[0:2]**2))
        if dictionary_SI['radius_2']<dist<dictionary_SI['radius_1']:
            b= dictionary_SI['rho_2']
        if dist<dictionary_SI['radius_2']:
            b= dictionary_SI['rho_1']
        if dist>dictionary_SI['radius_1']:
            b=0
        return b
#
#__________  radius_2
#rho_1     |
#          |
#          |
#          |
#          |         radius_1
#          |         ________________rho = 0
#          |        |
#          |        |
#          |        |
#          |________|rho_2 (note: rho_2 should be entered as a negative number if desired)
#
#



elif dictionary_SI['shape'] == 4:
    print "Monte Carlo Gaussian"
    def density(coords):
        b= np.exp(-(np.sum(coords[0:2]**2))/(dictionary_SI['radius_2']**2))
        return b

elif dictionary_SI['shape'] == 5:
    print "Chopped Cone Model - radius one is the radius closer to the -Z side of the box"
    def density(coords):
        #RHS is height of cone at the corrosponding z component.
        #LHS is radial distance of points
        if np.sqrt(np.sum(coords[0:2]**2)) < coords[2:3]*(dictionary_SI['radius_2']-dictionary_SI['radius_1'])/dictionary_SI['z_dim']+(dictionary_SI['radius_2']+dictionary_SI['radius_1'])/2:
            b = dictionary_SI['rho_1']
        else:
            b=0
        return b

elif dictionary_SI['shape'] == 6:
    print "Hexagonal Prism"
    def density(coords):
        coords = coords/dictionary_SI['radius_1']
        if coords[1]>0.5*3**0.5 or coords[1]<-0.5*3**0.5 or coords[1]+(coords[0]-1)*3**0.5>0 or coords[1]+(coords[0]+1)*3**0.5<0 or coords[1]-(coords[0]-1)*3**0.5<0 or coords[1]-(coords[0]+1)*3**0.5>0:
            b=0
        else:
            b=dictionary_SI['rho_1']
        return b

elif dictionary_SI['shape'] == 7:
    print "Rectangular prism"
    def density(coords):
        if coords[1:2]<dictionary_SI['radius_2']:
            return dictionary_SI['rho_1']
        else:
            return 0

elif dictionary_SI['shape'] == 8:
    print "String of bubbles"
    #distance between centres of spheres
    dist_centre = dictionary_SI['radius_2']
    def density(coords):
        #centre of nearest bubble:
        z_bubble = dist_centre/2 + dist_centre*np.floor(coords[2]/(dist_centre))
        #Disntance from the centre of this bubble
        dist = np.sqrt(np.sum((coords-[0,0,z_bubble])**2))
        if dist < dictionary_SI['radius_1']:
            b=dictionary_SI['rho_1']
        else:
            b=0
        return b

elif dictionary_SI['shape'] ==9:
    print "Chopped up cylinder"
    #gap between cylinders
    gap = dictionary_SI['radius_2']
    #number of gaps
    num_gaps = dictionary_SI['rho_2']
    #where does each gap begin?
    beginning_of_gap = [(random.random()-0.5)*dictionary_SI['z_dim'] for x in range(num_gaps)]   
    def density(coords):
        if np.sqrt(np.sum(coords[0:2]**2)) < dictionary_SI['radius_1'] and sum([1 for start_gap in beginning_of_gap if start_gap<coords[2]<start_gap+gap]) == 0:
            b=dictionary_SI['rho_1']
        else:
            b=0
        return b

elif dictionary_SI['shape'] == 10:
    print "custom radius depending on intensity"
    #the csv file MUST be one list. The z components are given by increments along the z axis.
    #Alternitively, get a list in anyway in the form:
    #custom = [1 2 3 4 5 6 7 8]
    z_dim = dictionary_SI['z_dim']

    custom = np.asarray(pylab.loadtxt(os.path.dirname(sys.argv[0])+"/custom.csv", delimiter=","))
    print custom
    def density(coords):
        
        low = np.floor(  (coords[2:3]+dictionary_SI['z_dim']/2)*(float(len(custom))-1.)/dictionary_SI['z_dim'])
        slope = (custom[int(low+1)] - custom[int(low)])*(float(len(custom))-1.)/dictionary_SI['z_dim']

        lin_interpolation = custom[int(low)]+ slope*(coords[2:3]-z_dim*(low/float(len(custom)-1)-0.5))
        #coords[2:3]*slope + custom[int(low)] -slope*dictionary_SI['z_dim']/(1.*len(custom)-1.)
        
        if lin_interpolation > np.sqrt(np.sum(coords[0:2]**2)):
            b = dictionary_SI['rho_1']
        else:
            b=0
        return b

elif dictionary_SI['shape'] == 11:
    print "double slit"
    
    def density(coords):
        if -dictionary_SI['radius_1']/2 <coords[0:1]<-dictionary_SI['radius_2']/2 or dictionary_SI['radius_2']/2 <coords[0:1]<dictionary_SI['radius_1']/2:
            return dictionary_SI['rho_1']
        else:
            return 0


