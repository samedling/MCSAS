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
   elif g.dictionary_SI['shape'] == 18:
      return d18doublecone_track(coords)
   elif g.dictionary_SI['shape'] == 19:
      return d19taperedcylinder(coords)
   elif g.dictionary_SI['shape'] == 20:
      return d20continuouscoreshell(coords)
   elif g.dictionary_SI['shape'] == 21:
      return d21coreshellcone(coords)
   elif g.dictionary_SI['shape'] == 22:
      return d22coreshellsmooth(coords) 
   elif g.dictionary_SI['shape'] == 23:
      g.dprint("Running density function when it should be importing points from file.")
      return d23importpoints(coords) 
      


#This is a dictionary containing all of the useful descriptions for the variables for each model.
#g.var_list=['radius_1','radius_2','z_dim','rho_1','rho_2','num','length_2']
g.model_parameters=[
   (0,'Analytic Model Only',('radius_1','radius_2','z_dim','rho_1','rho_2','num','length_2')),
   (1,'Sphere',("Radius (nm)","unused","unused","Density","unused","unused","unused")),
   (2,'Cylinder',("Radius (nm)","unused","Length (nm)","Density","unused","unused","unused")),
   (3,'Core Shell Cylinder',("Outer Radius (nm)","Inner Radius (nm)","Length (nm)","Core Density","Shell Density","unused","unused")),
   (4,'Gaussian Cylinder',("Radius (nm)","Std. Dev. (nm)","Length (nm)","Density","unused","unused","unused")),
   (5,'Chopped Cone',("Max Radius (nm)","Min Radius (nm)","Length (nm)","Density","unused","unused","unused")),
   (6,'Hexagonal Prism',("Side Length (nm)","unused","Length (nm)","Density","unused","unused","unused")),
   (7,'Rectangular Prism',("Long Side (nm)","Short Side (nm)","Length (nm)","Density","unused","unused","unused")),
   (8,'String of Bubbles',("Radius (nm)","Space Btwn Centers (nm)","Length (nm)","Density","unused","unused","unused")),
   (9,'Random Chopped Cylinder',("Radius (nm)","Gap Width (nm)","Length (nm)","Density","Number of Gaps","currently unused","unused")),
   (10,'CSV For Radius',('radius_1','radius_2','z_dim','rho_1','rho_2','num','length_2')),
   (11,'Double Slit',("Outside Distance (nm)","Inside Distance (nm)","Length (nm)","Density","unused","unused","unused")),
   (12,'N-gon Truncated Cone',("Radius for Large n (nm)","Radius for Small n (nm)","Length (nm)","Density","Number of Sides","currently unused","unused")),
   (13,'Sine Oscillation',("Origin to Peak","Origin to Trough","Length (nm)","Density","Number of Oscillations","currently unused","unused")),
   (14,'Double Cone',("End Radius (nm)","Central Radius (nm)","Length (nm)","Density","unused","unused","unused")),
   (15,'Elliptical Cylinder',("x Radius (nm)","y Radius (nm)","Length (nm)","Density","unused","unused","unused")),
   (16,'Aymm Hex Pyramid',("Orig. Side Length (nm)","Side Adjustment (nm)","Length (nm)","Density","unused","unused","unused")),
   (17,'Chopped Core Shell',("Outer Radius (nm)","Inner Radius (nm)","Length (nm)","Core Density","Shell Density","Number of Gaps","Gap Length (nm)")),
   (18,'Double Cone with Track',("End Cone Radius (nm)","Central Cone Radius (nm)","Total Length (nm)","Density","unused","unused","Track Radius (nm)")),
   (19,'Tapered Cylinder',('Radius (nm)','unused','Total Length (nm)','Density','unused','Taper Both Ends?','Cone Length (nm)')),
   (20,'Continuous Core Shell Cylinder',("Outer Radius (nm)","Inner Radius (nm)","Length (nm)","Extreme Core Density","Shell Density","unused","unused")),
   (21,'Core Shell Cone',("Outer Radius (nm)","Inner Radius (nm)","Length (nm)","Core Density","Shell Density","unused","Outer Radius - Far End")),
   (22,'Smooth Core Shell Cylinder',('Outer Radius (nm)','Core Radius (nm)','Length (nm)','Core Density','Shell Density','Smoothness (Low=smoother)','unused')),
   (23,'Load from "import.csv"', ('unused','unused','unused','unused','unused','unused','unused')),
]

#Template: fill in number, Name of Model, and useful descriptor or 'unused' in place of each variable name.
#(n,'Name of Model',('radius_1','radius_2','z_dim','rho_1','rho_2','num','length_2'))

g.MC_num_and_name = np.asarray([[g.model_parameters[i][1],i] for i in range(len(g.model_parameters))])


#This is the list of all the Monte Carlo Models that you choose from.
#g.MC_num_and_name = np.array([["Analytic Model Only",0],
#                        ["Sphere",1],
#                        ["Cylinder",2],
#                        ["Core shell cylinder",3],
#                        ["Gaussian",4],
#                        ["Cone Model",5],
#                        ["Hexagonal Prism",6],
#                        ["Rectangular Prism",7],
#                        ["String of Bubbles",8],
#                        ["Chopped up Cylinder",9],
#                        ["Custom CSV Defining the Radius",10],
#                        ["Double Slit",11],
#                        ["N-gon Truncated Cone",12],
#                        ["Sine Shaped Oscillation",13],
#                        ["Double Cone",14],
#                        ["Elliptical Cylinder",15],
#                        ["Asym Hex Pyramid",16],
#                        ["Chopped Core Shell",17]
#                        ])
g.MC_num_and_name_dict = {x[0]:x[1] for x in g.MC_num_and_name} #This is needed, so that when an option is chosen, we can find the shape number.


#def density_vector(all_coords):
#   '''Vectorizes any density function.'''
#   all_densities = np.empty(all_coords.shape[0])
#   for i in range(all_coords.shape[0]):
#      all_densities[i] = density(all_coords[i,:])
#   return all_densities


def d1sphere(coords):
    g.dictionary_SI['z_dim']=g.dictionary_SI['radius_1']*2.
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
   piece_l = (length-n_gaps*gap_l)/(n_gaps+1)
   g.dprint('Chopped Core Shell')
   print('Piece length is {0} nm.'.format(piece_l*10**9))

   if n_gaps%2:   #odd number of gaps
      #g.dprint('Odd number of gaps.')
      return [rho2 if r2<np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-gap_l/2)%(piece_l+gap_l) < piece_l else rho1 if np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-gap_l/2)%(piece_l+gap_l) < piece_l else 0 for i in range(coords.shape[0])]
   else:          #even number of gaps
      #g.dprint('Even number of gaps.')
      return [rho2 if r2<np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-piece_l/2)%(piece_l+gap_l) > gap_l else rho1 if np.sqrt(np.sum(coords[i,0:2]**2))<r1 and (np.abs(coords[i,2])-piece_l/2)%(piece_l+gap_l) > gap_l else 0 for i in range(coords.shape[0])]

def d18doublecone_track(coords):
    r1=g.dictionary_SI['radius_1']
    r2=g.dictionary_SI['radius_2']
    z_dim=g.dictionary_SI['z_dim']
    if g.dictionary_SI['radius_2'] < 0:
        print('Individual cone height is {0} nm.'.format(10**9*r1/(r1-r2)*z_dim/2))
        print('Space between cones is is {0} nm.'.format(10**9*-r2/(r1-r2)*z_dim))
    return [g.dictionary_SI['rho_1'] if ((np.sqrt(np.sum(coords[i,0:2]**2)) < np.abs(coords[i,2])*(g.dictionary_SI['radius_1']-g.dictionary_SI['radius_2'])/(g.dictionary_SI['z_dim']/2.)+g.dictionary_SI['radius_2']) or (np.sqrt(np.sum(coords[i,0:2]**2)) < g.dictionary_SI['length_2']) ) else 0 for i in range(coords.shape[0])]

def d19taperedcylinder(coords):
    radius = g.dictionary_SI['radius_1']
    rho = g.dictionary_SI['rho_1']
    length = g.dictionary_SI['z_dim']
    length_2 = g.dictionary_SI['length_2']
    if g.dictionary_SI['num']: #taper both ends
       return [rho if ((np.sqrt(np.sum(coords[i,0:2]**2)) < radius) and (np.sqrt(np.sum(coords[i,0:2]**2)) < -radius*np.abs(coords[i,2])/length_2+0.5*radius*length/length_2)) else 0 for i in range(coords.shape[0])]
    else: #taper only one end
       return [rho if ((np.sqrt(np.sum(coords[i,0:2]**2)) < radius) and (np.sqrt(np.sum(coords[i,0:2]**2)) < -radius*coords[i,2]/length_2+0.5*radius*length/length_2)) else 0 for i in range(coords.shape[0])]

def d20continuouscoreshell(coords):
    return [g.dictionary_SI['rho_2'] if g.dictionary_SI['radius_2']<np.sqrt(np.sum(coords[i,0:2]**2))<g.dictionary_SI['radius_1'] else
    (g.dictionary_SI['rho_2']-g.dictionary_SI['rho_1']/g.dictionary['radius_2'])*np.sqrt(np.sum(coords[i,0:2]**2))+g.dictionary_SI['rho_1']
    if np.sqrt(np.sum(coords[i,0:2]**2))<g.dictionary_SI['radius_2'] else 0 for i in range(coords.shape[0])]

def d21coreshellcone(coords):
    return [g.dictionary_SI['rho_1'] if np.sqrt(np.sum(coords[i,0:2]**2)) < coords[i,2:3]*(g.dictionary_SI['length_2']-g.dictionary_SI['radius_1'])*g.dictionary_SI['radius_2']/(g.dictionary_SI['z_dim']*g.dictionary_SI['radius_1'])+g.dictionary_SI['radius_2']+(g.dictionary_SI['length_2']-g.dictionary_SI['radius_1'])*g.dictionary_SI['radius_2']/(2*g.dictionary_SI['radius_1']) else g.dictionary_SI['rho_2'] if np.sqrt(np.sum(coords[i,0:2]**2)) < coords[i,2:3]*(g.dictionary_SI['length_2']-g.dictionary_SI['radius_1'])/g.dictionary_SI['z_dim']+(g.dictionary_SI['length_2']+g.dictionary_SI['radius_1'])/2 else 0 for i in range(coords.shape[0])]

def d22coreshellsmooth(coords):
    g.dictionary_SI['x_dim']=1.5*g.dictionary_SI['x_dim']
    g.dictionary_SI['y_dim']=1.5*g.dictionary_SI['y_dim']
    return [((g.dictionary_SI['rho_1']-g.dictionary_SI['rho_2'])*np.exp(-(np.sqrt(np.sum(coords[i,0:2]**2))/g.dictionary_SI['radius_2'])**g.dictionary_SI['num'])+ g.dictionary_SI['rho_2']*np.exp(-(np.sqrt(np.sum(coords[i,0:2]**2))/(g.dictionary_SI['radius_1']))**g.dictionary_SI['num'])) if np.sqrt(np.sum(coords[i,0:2]**2))<1.5*g.dictionary_SI['radius_1'] else 0 for i in range(coords.shape[0])]    

def d23importpoints(coords):
   None #In Monte_Carlo_Functions.py, we import points directly.

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






