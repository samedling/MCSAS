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

def compare(points):
   print('\n')
   #points = Points_For_Calculation()

   #Fast Way (check!)
   print('New Way (to check):')
   calc_int = Calculate_Intensity(points)
   #print(calc_int.sum())
   norm_calc_int = normalize(calc_int)
   #print norm_calc_int

   #Old Way
   print('Old Way (Infinite Coherence Length):')
   det_int = Detector_Intensity(points)  #runs with infinite coherence length
   #print(det_int.sum())
   norm_det_int = normalize(det_int)
   #print norm_det_int
   
   #Long Way
   print('Long Way:')
   long_int = long_way(points)
   #print('finished')
   #print(long_int.sum())
   norm_long_int = normalize(long_int)
   #print norm_long_int

   print('\n')
   old_diff = norm_long_int - norm_det_int
   new_diff = norm_long_int - norm_calc_int
   print('Old calculation differs by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(old_diff)),np.sum(old_diff**2)))
   print('New calculation differs actual by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(new_diff)),np.sum(new_diff**2)))

   #print('Plotting New Way, Old Way, and Long Way.')
   #Fit_plot(norm_calc_int,norm_det_int,norm_long_int)

   print('Plotting Long Way, New Way, and difference.')
   Fit_plot(norm_long_int,norm_calc_int,np.abs(norm_long_int-norm_calc_int))

   #print('Plotting Long Way, Old Way, and difference.')
   #Fit_plot(norm_long_int,norm_det_int,np.abs(norm_long_int-norm_det_int))

   print('All Done!')




if __name__ == '__main__':
   print('\n')
   g.dictionary['pixels'] = "50 50"
   g.dictionary['z_dim'] = 500
   g.dictionary['z_scale'] = 50
   #g.dictionary['d_lambda'] = 0     #infinite coherence length
   #g.dictionary['d_lambda'] = 2e-4  #correct coherence length (565nm)
   #g.dictionary['d_lambda'] = 5e-4  #coherence length = 225nm
   g.dictionary['d_lambda'] = 10e-4  #coherence length = 110nm
   #g.dictionary['d_lambda'] = 20e-4  #coherence length = 55nm


   make_SI_dict()
   
   points = Points_For_Calculation()
   compare(points)

#normalize & then take the difference?
