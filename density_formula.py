#THIS IS THE MOST IMPORTANT LINE!!!!!!!!!!!!!!!!
#When making points, x<x_dim/2, y<y_dim/2. When defining a function, by default, the x dimensions and y dimensions are equal to double radius_1, i.e. 2*radius_1 = x_dim = y_dim.
#Make sure this holds true for your function.
#you can always override x_dim and y_dim by defining them in your program (e.g. g.dictionary_SI['x_dim']=5*10**-9)
#This will NOT!! work for a sequence though.

#A reminder: the camera is a long way away in the z direction
#The function is run once for each point, which is entered as a numpy array, i.e. it is in the form: np.array([x,y,z])
#Make sure you always use g.dictionary_SI to get variables.

import pylab, sys, os, random
import numpy as np

import global_vars as g

def density_vector(all_coords):
   '''Vectorizes any density function.'''
   all_densities = np.empty(all_coords.shape[0])
   for i in range(all_coords.shape[0]):
      all_densities[i] = density(all_coords[i,:])
   return all_densities

if g.dictionary_SI['shape'] == 0:
    print "Analytic Model Only"

elif g.dictionary_SI['shape'] == 1:
    print "Monte Carlo Sphere"
    def density(coords):
        dist = np.sqrt(np.sum(coords**2))
        if dist < g.dictionary_SI['radius_1']:
            b=g.dictionary['rho_1']
        else:
            b=0
        return b

elif g.dictionary_SI['shape'] == 2:
    print "Monte Carlo Cylinder"
    def density(coords):
        dist = np.sqrt(np.sum(coords[0:2]**2))
        if dist < g.dictionary_SI['radius_1']:
            b=g.dictionary_SI['rho_1']
        else:
            b=0
        return b

elif g.dictionary_SI['shape'] == 3:
    print "Monte Carlo core shell cylinder"
    print "radius_2 must be less than radius_1"
    print "rho_1 is the density in the core (r<radius_2)"
    print "rho_2 is the density in the shell (radius_2<r<radius_1)"
    def density(coords):
        dist = np.sqrt(np.sum(coords[0:2]**2))
        if g.dictionary_SI['radius_2']<dist<g.dictionary_SI['radius_1']:
            b= g.dictionary_SI['rho_2']
        elif dist<g.dictionary_SI['radius_2']:
            b= g.dictionary_SI['rho_1']
        elif dist>g.dictionary_SI['radius_1']:
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

elif g.dictionary_SI['shape'] == 4:
    print "Monte Carlo Gaussian"
    def density(coords):
        return np.exp(-(np.sum(coords[0:2]**2))/(g.dictionary_SI['radius_2']**2))

elif g.dictionary_SI['shape'] == 5:
    print "Chopped Cone Model - radius one is the radius closer to the -Z side of the box"
    def density(coords):
        #RHS is height of cone at the corrosponding z component.
        #LHS is radial distance of points
        if np.sqrt(np.sum(coords[0:2]**2)) < coords[2:3]*(g.dictionary_SI['radius_2']-g.dictionary_SI['radius_1'])/g.dictionary_SI['z_dim']+(g.dictionary_SI['radius_2']+g.dictionary_SI['radius_1'])/2:
            b = g.dictionary_SI['rho_1']
        else:
            b=0
        return b

elif g.dictionary_SI['shape'] == 6:
    print "Hexagonal Prism"
    def density(coords):
        coords = coords/g.dictionary_SI['radius_1']
        if coords[1]>0.5*3**0.5 or coords[1]<-0.5*3**0.5 or coords[1]+(coords[0]-1)*3**0.5>0 or coords[1]+(coords[0]+1)*3**0.5<0 or coords[1]-(coords[0]-1)*3**0.5<0 or coords[1]-(coords[0]+1)*3**0.5>0:
            b=0
        else:
            b=g.dictionary_SI['rho_1']
        return b

elif g.dictionary_SI['shape'] == 7:
    print "Rectangular prism"
    def density(coords):
        if coords[1:2]<g.dictionary_SI['radius_2']:
            return g.dictionary_SI['rho_1']
        else:
            return 0

elif g.dictionary_SI['shape'] == 8:
    print "String of bubbles"
    #distance between centres of spheres
    dist_centre = g.dictionary_SI['radius_2']
    def density(coords):
        #centre of nearest bubble:
        z_bubble = dist_centre/2 + dist_centre*np.floor(coords[2]/(dist_centre))
        #Disntance from the centre of this bubble
        dist = np.sqrt(np.sum((coords-[0,0,z_bubble])**2))
        if dist < g.dictionary_SI['radius_1']:
            b=g.dictionary_SI['rho_1']
        else:
            b=0
        return b

elif g.dictionary_SI['shape'] ==9:
    print "Chopped up cylinder"
    #gap between cylinders
    gap = g.dictionary_SI['radius_2']
    #number of gaps
    num_gaps = g.dictionary_SI['rho_2']
    #where does each gap begin?
    beginning_of_gap = [(random.random()-0.5)*g.dictionary_SI['z_dim'] for x in range(num_gaps)]   
    def density(coords):
        if np.sqrt(np.sum(coords[0:2]**2)) < g.dictionary_SI['radius_1'] and sum([1 for start_gap in beginning_of_gap if start_gap<coords[2]<start_gap+gap]) == 0:
            b=g.dictionary_SI['rho_1']
        else:
            b=0
        return b

elif g.dictionary_SI['shape'] == 10:
    print "custom radius depending on intensity"
    #the csv file MUST be one list. The z components are given by increments along the z axis.
    #Alternitively, get a list in anyway in the form:
    #custom = [1 2 3 4 5 6 7 8]
    z_dim = g.dictionary_SI['z_dim']

    custom = np.asarray(pylab.loadtxt(os.path.dirname(sys.argv[0])+"/custom.csv", delimiter=","))
    print custom
    def density(coords):
        
        low = np.floor(  (coords[2:3]+g.dictionary_SI['z_dim']/2)*(float(len(custom))-1.)/g.dictionary_SI['z_dim'])
        slope = (custom[int(low+1)] - custom[int(low)])*(float(len(custom))-1.)/g.dictionary_SI['z_dim']

        lin_interpolation = custom[int(low)]+ slope*(coords[2:3]-z_dim*(low/float(len(custom)-1)-0.5))
        #coords[2:3]*slope + custom[int(low)] -slope*g.dictionary_SI['z_dim']/(1.*len(custom)-1.)
        
        if lin_interpolation > np.sqrt(np.sum(coords[0:2]**2)):
            b = g.dictionary_SI['rho_1']
        else:
            b=0
        return b

elif g.dictionary_SI['shape'] == 11:
    print "double slit"
    def density(coords):
        if -g.dictionary_SI['radius_1']/2 <coords[0:1]<-g.dictionary_SI['radius_2']/2 or g.dictionary_SI['radius_2']/2 <coords[0:1]<g.dictionary_SI['radius_1']/2:
            return g.dictionary_SI['rho_1']
        else:
            return 0

elif g.dictionary_SI['shape']==12:
    print "N-gon Truncated Cone"
    #rho_2 is the number of sides
    def density(coords):
        angle_number = np.float(np.floor(np.angle(coords[0:1]+1j*coords[1:2])*g.dictionary_SI['rho_2']/(2*3.14159265)))
        zrad = coords[2:3]*(g.dictionary_SI['radius_2']-g.dictionary_SI['radius_1'])/g.dictionary_SI['z_dim']+(g.dictionary_SI['radius_2']+g.dictionary_SI['radius_1'])/2.
    
        slope = (np.sin((angle_number+1)*2.*3.14159265/float(g.dictionary_SI['rho_2']))-np.sin((angle_number)*2.*3.14159265/float(g.dictionary_SI['rho_2'])))/(np.cos((angle_number+1)*2.*3.14159265/float(g.dictionary_SI['rho_2']))-np.cos((angle_number)*2.*3.14159265/float(g.dictionary_SI['rho_2'])))

        x_intercept = zrad*(slope*np.cos(angle_number*2*3.14159265/float(g.dictionary_SI['rho_2']))-np.sin(angle_number*2*3.14159265/float(g.dictionary_SI['rho_2'])))/(coords[1:2]/coords[0:1]-slope)
        y_intercept = x_intercept*coords[1:2]/coords[0:1]
        if x_intercept**2+y_intercept**2 > coords[0:1]**2+coords[1:2]**2:
            return g.dictionary_SI['rho_1']
        else:
            return 0

elif g.dictionary_SI['shape'] == 13:
    print "Sine-shaped-oscillation"
    def density(coords):
        if np.sqrt(np.sum(coords[0:2]**2)) < (g.dictionary_SI['radius_1']+g.dictionary_SI['radius_2'])/2 + (g.dictionary_SI['radius_1']-g.dictionary_SI['radius_2'])/2*np.sin(coords[2:3]*g.dictionary_SI['rho_2']*2*np.pi/g.dictionary_SI['z_dim']):
            b = g.dictionary_SI['rho_1']
        else:
            b=0
        return b            
           
elif g.dictionary_SI['shape'] == 14:
    print 'Double Cone'
    def density(coords):
        if np.sqrt(np.sum(coords[0:2]**2)) < abs(coords[2:3])*(g.dictionary_SI['radius_1']-g.dictionary_SI['radius_2'])/(g.dictionary_SI['z_dim']/2.)+g.dictionary_SI['radius_2']:
            b = g.dictionary_SI['rho_1']
        else:
            b=0
        return b  
 
elif g.dictionary_SI['shape'] == 15:
    print("Eliptical Cylinder")
    print("radius 1 is the x-axis; radius 2 is the y-axis.")
    def density(coords):
        dist = np.sqrt(coords[0]**2/g.dictionary_SI['radius_1']**2 + coords[1]**2/g.dictionary_SI['radius_2']**2)
        if dist < 1:
            b=g.dictionary_SI['rho_1']
        else:
            b=0
        return b



