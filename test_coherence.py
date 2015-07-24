#!/usr/bin/python

import numpy as np
from newgui import *



def compare(points,speed='Fortran'):
   print('{0}: Starting calculations...'.format(time.strftime("%X")))
   #points = Points_For_Calculation()

   #Fast Way (check!)
   print('New Way (to check):')
   calc_int = Calculate_Intensity(points,coherence_dup=1,coherence_taper=False)
   #print(calc_int.sum())
   calc_int = normalize(calc_int)
   #print norm_calc_int

   #Old Way
   print('Old Way (Infinite Coherence Length):')
   det_int = Detector_Intensity(points)  #runs with infinite coherence length
   #print(det_int.sum())
   det_int = normalize(det_int)
   #print norm_det_int
   
   #Long Way
   print('Long Way:')
   print('{0}: Starting long calculations using {1}...'.format(time.strftime("%X"),speed))
   long_int = Accurate_Intensity(points,speed)
   #print('finished')
   #print(long_int.sum())
   long_int = normalize(long_int)
   #print norm_long_int

   print('\n')
   old_diff = long_int - det_int
   new_diff = long_int - calc_int
   print('Old calculation differs by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(old_diff)),np.sum(old_diff**2)))
   print('New calculation differs actual by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(new_diff)),np.sum(new_diff**2)))


   #plot_all()
   #plot_new()
   #plot_old()

   return calc_int,det_int,long_int

def print_diff(test,correct,name=''):
   diff = correct - test
   print('{2} Calculation differs by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(diff)),np.sum(diff**2),name))

def plot_all():
   print('Plotting New Way, Old Way, and Long Way.')
   Fit_plot(g.calc_int,g.det_int,g.long_int)

def plot_new():
   print('Plotting Long Way, New Way, and difference.')
   Fit_plot(g.long_int,g.calc_int,np.abs(g.long_int-g.calc_int))

def plot_old():
   print('Plotting Long Way, Old Way, and difference.')
   Fit_plot(g.long_int,g.det_int,np.abs(g.long_int-g.det_int))


def average_test(times=5):
   print('{0}: Starting Calculation 1.'.format(time.strftime("%X")))
   g.points = Points_For_Calculation()
   g.calc_int,g.det_int,g.long_int = compare(g.points,'Fortran')
   for i in range(times-1):
      print('\n')
      print('{0}: Starting Calculation {1}.'.format(time.strftime("%X"),i+2))
      g.points = Points_For_Calculation()
      calc_int,det_int,long_int = compare(g.points,'Fortran')
      g.calc_int += calc_int
      g.det_int += det_int
      g.long_int += long_int
   g.calc_int = normalize(g.calc_int)
   g.det_int = normalize(g.det_int)
   g.long_int = normalize(g.long_int)
   print('\n')
   print_diff(g.det_int,g.long_int,name='Old')
   print_diff(g.calc_int,g.long_int,name='New')



if __name__ == '__main__':
   g.verbose = False
   np.random.seed([2015])
   print('\n')
   g.dictionary['pixels'] = "50 50"
   g.dictionary['z_dim'] = 250
   g.dictionary['z_scale'] = 50
   #g.dictionary['d_lambda'] = 0     #infinite coherence length
   #g.dictionary['d_lambda'] = 2e-4  #correct coherence length (565nm)
   g.dictionary['d_lambda'] = 5e-4  #coherence length = 225nm
   #g.dictionary['d_lambda'] = 10e-4  #coherence length = 110nm
   #g.dictionary['d_lambda'] = 20e-4  #coherence length = 55nm

   make_SI_dict()
   
   try:
      if not len(g.points):
         g.points = Points_For_Calculation()
         g.accurate = Accurate_Intensity(g.points)
   except AttributeError:
      g.points = Points_For_Calculation()
      g.accurate = normalize(Accurate_Intensity(g.points))
   dup = (1,2,5,10)
   taper = (0,1)
   for i in range(len(dup)*len(taper)):
      test = Calculate_Intensity(g.points,coherence_dup=dup[i/2],coherence_taper=taper[i%2])
      #diff = g.accurate - normalize(test)
      diff = normalize(g.accurate) - normalize(test)
      print('Calculation ({2},{3}) differs by {0:6.4f} (least sq: {1:8.6f}).'.format(np.sum(np.abs(diff)),np.sum(diff**2),dup[i/2],taper[i%2]))


   print('{0}: All Done!'.format(time.strftime("%X")))

   #plot_all()

#normalize & then take the difference?
