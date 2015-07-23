#!/usr/bin/python

import numpy as np
from newgui import *



#Long Way
def long_way(Points,speed='fortran'):
   symmetric = g.dictionary_SI['symmetric']
   QSize = g.dictionary_SI['QSize']
   x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
   EHC = g.dictionary_SI['EHC']
   mask = np.ones((y_pixels,x_pixels))
   coherence_length = 2*np.pi/(g.dictionary_SI['EHC']*g.dictionary['d_lambda'])
   #print('Object length is {0}'.format(Points[-1,2]-Points[0,2])) #wrong; points not sorted
   #print('Coherence length is {0}'.format(coherence_length))
   if speed[0] in ('f','F'):  #fortran
      #print('Calling Fortran...')
      return fastmath.sumint.coherent_long(QSize,EHC,coherence_length,mask,Points.T)
   elif speed[0] in ('o','O'):   #OpenCL
      #print('Calling OpenCL...')
      return g.opencl_sumint.sumint(QSize,EHC,x_pixels,y_pixels,Points,symmetric,coherence_length=coherence_length)
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

def compare(points,speed='Fortran'):
   print('\n')
   print('{0}: Starting calculations...'.format(time.strftime("%X"),speed))
   #points = Points_For_Calculation()

   #Fast Way (check!)
   print('New Way (to check):')
   calc_int = Calculate_Intensity(points,coherence_dup=1,coherence_taper=False)
   #print(calc_int.sum())
   g.calc_int = normalize(calc_int)
   #print norm_calc_int

   #Old Way
   print('Old Way (Infinite Coherence Length):')
   det_int = Detector_Intensity(points)  #runs with infinite coherence length
   #print(det_int.sum())
   g.det_int = normalize(det_int)
   #print norm_det_int
   
   #Long Way
   print('Long Way:')
   print('{0}: Starting long calculations using {1}...'.format(time.strftime("%X"),speed))
   long_int = long_way(points,speed)
   #print('finished')
   #print(long_int.sum())
   g.long_int = normalize(long_int)
   #print norm_long_int

   print('\n')
   old_diff = g.long_int - g.det_int
   new_diff = g.long_int - g.calc_int
   print('Old calculation differs by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(old_diff)),np.sum(old_diff**2)))
   print('New calculation differs actual by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(new_diff)),np.sum(new_diff**2)))

   print('{0}: All Done!'.format(time.strftime("%X")))

   plot_all()
   #plot_new()
   #plot_old()

def plot_all():
   print('Plotting New Way, Old Way, and Long Way.')
   Fit_plot(g.calc_int,g.det_int,g.long_int)

def plot_new():
   print('Plotting Long Way, New Way, and difference.')
   Fit_plot(g.long_int,g.calc_int,np.abs(g.long_int-g.calc_int))

def plot_old():
   print('Plotting Long Way, Old Way, and difference.')
   Fit_plot(g.long_int,g.det_int,np.abs(g.long_int-g.det_int))




if __name__ == '__main__':
   g.verbose = False
   print('\n')
   g.dictionary['pixels'] = "20 20"
   g.dictionary['z_dim'] = 500
   g.dictionary['z_scale'] = 50
   #g.dictionary['d_lambda'] = 0     #infinite coherence length
   #g.dictionary['d_lambda'] = 2e-4  #correct coherence length (565nm)
   #g.dictionary['d_lambda'] = 5e-4  #coherence length = 225nm
   g.dictionary['d_lambda'] = 10e-4  #coherence length = 110nm
   #g.dictionary['d_lambda'] = 20e-4  #coherence length = 55nm


   make_SI_dict()
   
   g.points = Points_For_Calculation()
   compare(g.points,'OpenCL')

#normalize & then take the difference?
