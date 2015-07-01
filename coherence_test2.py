#!/usr/bin/python

import numpy as np
from newgui import *
g.dictionary['pixels'] = "20 20"
g.dictionary['z_dim'] = 100
g.dictionary['z_scale'] = 10
g.dictionary['d_lambda'] = 5e-4   #0 or 2e-4
make_SI_dict()
points = Points_For_Calculation()

#Old Way
det_int = Detector_Intensity(points)  #runs with infinite coherence length
print(det_int.sum())
norm_det_int = normalize(det_int)


#Long Way
def long_way(Points):
   symmetric = g.dictionary_SI['symmetric']
   QSize = g.dictionary_SI['QSize']
   x_pixels,y_pixels = [int(i) for i in g.dictionary_SI['pixels'].split()]
   EHC = g.dictionary_SI['EHC']
   mask = np.ones((y_pixels,x_pixels))
   coherence_length = 2*np.pi/(g.dictionary_SI['EHC']*g.dictionary['d_lambda'])
   return fastmath.sumint.coherent_long(QSize,EHC,coherence_length,mask,Points.T)

long_int = long_way(points)
print('finished')
print(long_int.sum())
norm_long_int = normalize(long_int)

#Fast Way (check!)
#calc_int = Calculate_Intensity(points)
#print(calc_int.sum())
#norm_calc_int = normalize(calc_int)


#normalize & then take the difference?
