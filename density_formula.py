#THIS IS THE MOST IMPORTANT LINE!!!!!!!!!!!!!!!!
#When making points, x<x_dim/2, y<y_dim/2. When defining a function, by default, the x dimensions and y dimensions are equal to double radius_1, i.e. 2*radius_1 = x_dim = y_dim.
#Make sure this holds true for your function.
#you can always override x_dim and y_dim by defining them in your program (e.g. g.dictionary_SI['x_dim']=5*10**-9)
#This will NOT!! work for a sequence though.

#A reminder: the camera is a long way away in the z direction
#The function is run once for each point, which is entered as a numpy array, i.e. it is in the form: np.array([x,y,z])
#Make sure you always use g.dictionary_SI to get variables.

#TODO: vectorize rest and move original code next to list comprehension versions as comments

import pylab, sys, os, random
import numpy as np

import global_vars as g


def density(coords):
   if g.dictionary_SI['shape'] == 1:
      return d1sphere(coords)
   elif g.dictionary_SI['shape'] == 2:
      return d2cylinder(coords)
   elif g.dictionary_SI['shape'] == 3:
      return d3coreshell(coords)
   elif g.dictionary_SI['shape'] == 4:
      return d4gaussian(coords)
   elif g.dictionary_SI['shape'] == 5:
      return d5choppedcone(coords)
   elif g.dictionary_SI['shape'] == 6:
      return d6hexprism(coords)
   elif g.dictionary_SI['shape'] == 7:
      return d7rectprism(coords)
   elif g.dictionary_SI['shape'] == 8:
      return d8bubbles(coords)
   elif g.dictionary_SI['shape'] == 9:
      return d9choppedcyl(coords)
   elif g.dictionary_SI['shape'] == 10:
      return d10custom(coords)
   elif g.dictionary_SI['shape'] == 11:
      return d11doubleslit(coords)
   elif g.dictionary_SI['shape'] == 12:
      return d12ngon(coords)
   elif g.dictionary_SI['shape'] == 13:
      return d13sine(coords)
   elif g.dictionary_SI['shape'] == 14:
      return d14doublecone(coords)
   elif g.dictionary_SI['shape'] == 15:
      return d15ellipticalcylinder(coords)
   elif g.dictionary_SI['shape'] == 16:
      return d16asymmhexpyr(coords)
   elif g.dictionary_SI['shape'] == 17:
      return d17choppedcoreshell(coords)



#def density_vector(all_coords):
#   '''Vectorizes any density function.'''
#   all_densities = np.empty(all_coords.shape[0])
#   for i in range(all_coords.shape[0]):
#      all_densities[i] = density(all_coords[i,:])
#   return all_densities


def d1sphere(coords):
    return [g.dictionary_SI['rho_1'] if np.sqrt(np.sum(coords[i,:]**2)) < g.dictionary_SI['radius_1'] else 0 for i in range(coords.shape[0])]

def d2cylinder(coords):
    return [g.dictionary_SI['rho_1'] if np.sqrt(np.sum(coords[i,0:2]**2)) < g.dictionary_SI['radius_1'] else 0 for i in range(coords.shape[0])]

def d3coreshell(coords):
    return [g.dictionary_SI['rho_2'] if g.dictionary_SI['radius_2']<np.sqrt(np.sum(coords[i,0:2]**2))<g.dictionary_SI['radius_1'] else g.dictionary_SI['rho_1'] if np.sqrt(np.sum(coords[i,0:2]**2))<g.dictionary_SI['radius_2'] else 0 for i in range(coords.shape[0])]

def d4gaussian(coords):
    return [np.exp(-(np.sum(coords[i,0:2]**2))/(g.dictionary_SI['radius_2']**2)) for i in range(coords.shape[0])]

def d5choppedcone(coords):
    return [g.dictionary_SI['rho_1'] if np.sqrt(np.sum(coords[i,0:2]**2)) < coords[i,2:3]*(g.dictionary_SI['radius_2']-g.dictionary_SI['radius_1'])/g.dictionary_SI['z_dim']+(g.dictionary_SI['radius_2']+g.dictionary_SI['radius_1'])/2 else 0 for i in range(coords.shape[0])]

def d6hexprism(coords):
    #r1 = g.dictionary_SI['radius_1']
    coords = coords/g.dictionary_SI['radius_1']
    return [0 if (coords[i,1]>np.sqrt(3)/2 or coords[i,1]<-np.sqrt(3)/2 or coords[i,1]+(coords[i,0]-1)*np.sqrt(3)>0 or coords[i,1]+(coords[i,0]+1)*np.sqrt(3)<0 or coords[i,1]-(coords[i,0]-1)*np.sqrt(3)<0 or coords[i,1]-(coords[i,0]+1)*np.sqrt(3)>0) else g.dictionary_SI['rho_1'] for i in range(coords.shape[0])]

def d7rectprism(coords):
    return [g.dictionary_SI['rho_1'] if coords[i,1]<g.dictionary_SI['radius_2'] else 0 for i in range(coords.shape[0])]

def d8bubbles(coords):
    dist_centre = g.dictionary_SI['radius_2']
    return [g.dictionary_SI['rho_1'] if np.sqrt(np.sum((coords[i,:]-[0,0,dist_centre/2+dist_centre*np.floor(coords[i,2]/(dist_centre))])**2)) < g.dictionary_SI['radius_1'] else 0 for i in range(coords.shape[0])]

def d9choppedcyl(coords):
    #todo: improve algorithm?
    #TODO: random number generator prevents fitting from working!!
    gap = g.dictionary_SI['radius_2']
    num_gaps = g.dictionary_SI['rho_2']
    return [g.dictionary_SI['rho_1'] if (np.sqrt(np.sum(coords[0:2]**2)) < g.dictionary_SI['radius_1'] and np.sum([1 for start_gap in [(random.random()-0.5)*g.dictionary_SI['z_dim'] for x in range(num_gaps)] if start_gap<coords[2]<start_gap+gap]) == 0) else 0 for i in range(coords.shape[0])]

def d10custom(coords):
    z_dim = g.dictionary_SI['z_dim']
    custom = np.asarray(pylab.loadtxt(os.path.dirname(sys.argv[0])+"/custom.csv", delimiter=","))
    densities = np.empty(coords.shape[0])

    for i in range(coords.shape[0]):
        
        low = np.floor(  (coords[i,2]+g.dictionary_SI['z_dim']/2)*(float(len(custom))-1.)/g.dictionary_SI['z_dim'])
        slope = (custom[int(low+1)] - custom[int(low)])*(float(len(custom))-1.)/g.dictionary_SI['z_dim']

        lin_interpolation = custom[int(low)]+ slope*(coords[i,2]-z_dim*(low/float(len(custom)-1)-0.5))
        #coords[2:3]*slope + custom[int(low)] -slope*g.dictionary_SI['z_dim']/(1.*len(custom)-1.)
        
        if lin_interpolation > np.sqrt(np.sum(coords[i,0:2]**2)):
            densities[i] = g.dictionary['rho_1']
        else:
            densities[i] = 0
    return densities

def d11doubleslit(coords):
    return [g.dictionary_SI['rho_1'] if (-g.dictionary_SI['radius_1']/2 <coords[i,0]<-g.dictionary_SI['radius_2']/2 or g.dictionary_SI['radius_2']/2 <coords[i,0]<g.dictionary_SI['radius_1']/2) else 0 for i in range(coords.shape[0])]

def d12ngon(coords):
    densities = np.empty(coords.shape[0])
    for i in range(coords.shape[0]):
        angle_number = np.float(np.floor(np.angle(coords[i,0]+1j*coords[i,1])*g.dictionary_SI['rho_2']/(2*np.pi)))
        zrad = coords[i,2]*(g.dictionary_SI['radius_2']-g.dictionary_SI['radius_1'])/g.dictionary_SI['z_dim']+(g.dictionary_SI['radius_2']+g.dictionary_SI['radius_1'])/2.
        
        slope = (np.sin((angle_number+1)*2.*np.pi/float(g.dictionary_SI['rho_2']))-np.sin((angle_number)*2.*np.pi/float(g.dictionary_SI['rho_2'])))/(np.cos((angle_number+1)*2.*np.pi/float(g.dictionary_SI['rho_2']))-np.cos((angle_number)*2.*np.pi/float(g.dictionary_SI['rho_2'])))

        x_intercept = zrad*(slope*np.cos(angle_number*2*np.pi/float(g.dictionary_SI['rho_2']))-np.sin(angle_number*2*np.pi/float(g.dictionary_SI['rho_2'])))/(coords[i,1]/coords[i,0]-slope)
        y_intercept = x_intercept*coords[i,1]/coords[i,0]
        if x_intercept**2+y_intercept**2 > coords[i,0]**2+coords[i,1]**2:
            densities[i] = g.dictionary['rho_1']
        else:
            densities[i] = 0
    return densities

def d13sine(coords):
    return [g.dictionary_SI['rho_1'] if (np.sqrt(np.sum(coords[i,0:2]**2)) < (g.dictionary_SI['radius_1']+g.dictionary_SI['radius_2'])/2 + (g.dictionary_SI['radius_1']-g.dictionary_SI['radius_2'])/2*np.sin(coords[i,2]*g.dictionary_SI['rho_2']*2*np.pi/g.dictionary_SI['z_dim'])) else 0 for i in range(coords.shape[0])]

def d14doublecone(coords):
    return [g.dictionary_SI['rho_1'] if (np.sqrt(np.sum(coords[i,0:2]**2)) < np.abs(coords[i,2])*(g.dictionary_SI['radius_1']-g.dictionary_SI['radius_2'])/(g.dictionary_SI['z_dim']/2.)+g.dictionary_SI['radius_2']) else 0 for i in range(coords.shape[0])]

def d15ellipticalcylinder(coords):
    return [g.dictionary_SI['rho_1'] if coords[i,0]**2/g.dictionary_SI['radius_1']**2 + coords[i,1]**2/g.dictionary_SI['radius_2']**2 < 1 else 0 for i in range(coords.shape[0])]

def d16asymmhexpyr(coords):
    r1=g.dictionary_SI['radius_1']
    r2=g.dictionary_SI['radius_2']
    offset = r2/r1
    densities = np.empty(coords.shape[0])
    for i in range(coords.shape[0]):
        scale_by = 2*g.dictionary_SI['z_dim'] / (g.dictionary_SI['z_dim']-2*coords[i,2]) * max(1,(r1+r2)/r1)/r1      #first couple term term makes a single cone, second term scales width so even if enlarged it fits within random points box
        x=coords[i,0]*scale_by
        y=coords[i,1]*scale_by
        if np.abs(y) > np.sqrt(3)/2 or np.abs(y) > (-np.sqrt(3)*(np.abs(x)-offset)+np.sqrt(3)):
            densities[i] = 0
        else:
            densities[i] = g.dictionary['rho_1']
    return densities

def d17choppedcoreshell(coords):
   r1 = g.dictionary_SI['radius_1']
   r2 = g.dictionary_SI['radius_2']
   rho1 = g.dictionary_SI['rho_1']
   rho2 = g.dictionary_SI['rho_2']
   length = g.dictionary_SI['z_dim']
   n_gaps = g.dictionary_SI['num']
   gap_l = g.dictionary_SI['length_2']
   piece_l = length-n_gaps*gap_l
   if n_gaps%2:   #odd number of gaps
      return [rho2 if r2<np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-piece_l/2)%(piece_l+gap_l) > gap_l else rho1 if np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-piece_l/2)%(piece_l+gap_l) > gap_l else 0 for i in range(coords.shape[0])]
   else:          #even number of gaps
      return [rho2 if r2<np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-gap_l/2)%(piece_l+gap_l) < piece_l else rho1 if np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-gap_l/2)%(piece_l+gap_l) < piece_l else 0 for i in range(coords.shape[0])]



#def density_slow_template(coords):
    #densities = np.empty(coords.shape[0])   #creates empty density array
    #for i in range(coords.shape[0]):        #iterates through all coordinates
        #if np.sqrt(np.sum(coords[i,:]**2)) < g.dictionary_SI['radius_1']:    #EDIT
            #densities[i] = g.dictionary['rho_1']                             #THESE
        #else:                                                                #FOUR
            #densities[i] = 0                                                 #LINES (and add additional elif statements if needed).
    #return densities                        #returns density array

#def density_fast_template(coords):
    #return [g.dictionary_SI['rho_1'] if <insert if statement here using coords[i,:]> else 0 for i in range(coords.shape[0])]






